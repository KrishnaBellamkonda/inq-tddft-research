/*
 * inqkit::SelfInteractionCorrection — remove the LDA self-interaction of ONE
 * tagged Kohn-Sham orbital (the projectile wavepacket) during real-time
 * propagation, wrapper-only (no inq/ or inq-study/ edit).
 *
 * Scheme (docs/plans/wp-self-interaction-correction.md, reviewed 2026-08-02):
 * between INQ propagator steps, the WP column of phi.matrix() receives a
 * multiplicative correction kick and is then projected back out of the
 * occupied manifold:
 *
 *     psi_wp  <-  N . Q . exp(+i dt_eff v_SIC) psi_wp
 *     Q  =  1 - sum_{j occupied} |psi_j><psi_j|
 *
 * with, per variant (Mode):
 *     hartree (SIC-H)   : v_SIC = v_H[n_wp]
 *     pz_run  (SIC-PZrun): v_SIC = v_H[n_wp] + v_xc^unpol[n_wp]
 *
 * n_wp = occ_wp |psi_wp|^2 (occ_wp asserted = 1). The xc self-term is
 * evaluated with the SAME spin treatment as the run (UNPOLARISED libxc
 * XC_LDA_X + XC_LDA_C_PZ, through INQ's own xc_term::evaluate_functional), so
 * that in a one-electron vacuum run the subtraction cancels the Hamiltonian's
 * own v_H + v_xc EXACTLY and the packet disperses freely — the closed-form
 * validation gate. Canonical PZ would use the POLARISED functional
 * (Perdew & Zunger, PRB 23, 5048 (1981), Eq. 34); in a spin-restricted run
 * that subtracts a term the Hamiltonian never contained (plan §0/D1).
 *
 * The projected kick is the ONE-SIDED Lagrange-multiplier form of variational
 * time-dependent SIC (Messud, Dinh, Reinhard, Suraud, PRL 101, 096404 (2008)):
 * it preserves WP-bath orthonormality EXACTLY (V|psi_j> = 0 kills the leak
 * rate d/dt<psi_wp|psi_j> = -i<psi_wp|V|psi_j>) but is NOT the full
 * variational scheme, so the corrected energy
 *     E_corr = E_KS - U[n_wp] (- E_xc^unpol[n_wp])
 * is exactly conserved only when no occupied state overlaps the packet
 * (vacuum). In the jellium its drift is a MEASURED diagnostic (plan §0/D2;
 * cf. Mundt et al., PRA 75, 050501(R) (2007) for why residuals of simplified
 * TD-SIC schemes must be watched for secular growth).
 *
 * Every step reports: u_self = U[n_wp], exc_self = E_xc^unpol[n_wp],
 * max_overlap_pre (post-kick, pre-projection bath overlap — the leak the
 * projection removed), norm_removed (fractional norm taken by Q), so the run
 * can write sic.csv and gate on it.
 *
 * Single k-point, single rank (same limits as inqkit::WavePacket).
 */
#pragma once

#include <inq/inq.hpp>

#include <cmath>
#include <stdexcept>
#include <string>

#ifdef __CUDACC__
#include <cuda_runtime.h>
#ifndef INQKIT_GPU_SYNC
#define INQKIT_GPU_SYNC() cudaDeviceSynchronize()
#endif
#else
#ifndef INQKIT_GPU_SYNC
#define INQKIT_GPU_SYNC() ((void)0)
#endif
#endif

namespace inqkit {

struct SICStepReport {
	double u_self = 0.0;           // U[n_wp] = 1/2 integral n_wp phi_wp   (Ha)
	double exc_self = 0.0;         // E_xc^unpol[n_wp]                     (Ha)
	double max_overlap_pre = 0.0;  // max_j |<psi_j|psi_wp>| after the kick, before Q
	double norm_removed = 0.0;     // fractional norm^2 removed by Q this step
	int    n_projected = 0;        // occupied states projected against
	bool   kicked = false;         // false in Mode::none or dt_eff == 0
};

class SelfInteractionCorrection {
public:
	enum class Mode { none, hartree, pz_run };

	static Mode mode_from_string(std::string const & s) {
		if(s == "none" || s == "")  return Mode::none;
		if(s == "h" || s == "hartree") return Mode::hartree;
		if(s == "pzrun" || s == "pz_run" || s == "pz") return Mode::pz_run;
		throw std::runtime_error("SelfInteractionCorrection: unknown mode '" + s +
		                         "' (want none|h|pzrun)");
	}
	static char const * mode_name(Mode m) {
		switch(m){ case Mode::none: return "none";
		           case Mode::hartree: return "sic-h";
		           default: return "sic-pzrun"; }
	}

	SelfInteractionCorrection(Mode mode, int wp_index, double occ_threshold = 1e-8)
		: mode_(mode), wp_(wp_index), occ_thr_(occ_threshold) {}

	Mode mode() const { return mode_; }

	// U[n_wp] and E_xc^unpol[n_wp] of the CURRENT state, no state change.
	// Valid in every mode (diagnostics for uncorrected runs too).
	template <class Electrons>
	SICStepReport measure(Electrons const & el) const {
		SICStepReport rep;
		auto n_wp = wp_density_(el);
		auto phi_wp = inq::solvers::poisson::solve(n_wp);
		rep.u_self = 0.5 * inq::operations::integral_product(n_wp, phi_wp);
		rep.exc_self = exc_only_(n_wp);
		return rep;
	}

	// The full correction step: build v_SIC from the CURRENT state, kick the WP
	// column by exp(+i dt_eff v_SIC), project out the occupied manifold,
	// renormalise. Returns the energies of the PRE-kick state (a real kick does
	// not change n_wp, so u_self/exc_self are kick-invariant) plus the
	// projection diagnostics. Mode::none measures and returns without touching
	// the state.
	template <class Electrons>
	SICStepReport apply(Electrons & el, double dt_eff) const {
		if(el.kpin().size() != 1)
			throw std::runtime_error("SelfInteractionCorrection: gamma-only runs only.");
		auto & phi = el.kpin()[0];
		if(phi.basis().comm().size() != 1 || phi.set_comm().size() != 1)
			throw std::runtime_error("SelfInteractionCorrection: single-rank only.");

		SICStepReport rep;

		// occ_wp = 1 is a precondition of the one-electron PZ logic.
		const double occ_wp = el.occupations()[0][wp_];
		if(mode_ != Mode::none && std::abs(occ_wp - 1.0) > 1e-6)
			throw std::runtime_error("SelfInteractionCorrection: WP occupation must be 1.");

		auto n_wp = wp_density_(el);
		auto v = inq::solvers::poisson::solve(n_wp);          // v_H[n_wp]
		rep.u_self = 0.5 * inq::operations::integral_product(n_wp, v);
		rep.exc_self = exc_only_(n_wp);

		if(mode_ == Mode::none || dt_eff == 0.0) return rep;

		if(mode_ == Mode::pz_run) add_vxc_unpol_(n_wp, v);    // v += v_xc^unpol[n_wp]

		const int n_pts = phi.basis().local_size();
		const double dV = phi.basis().volume_element();
		auto mat_ = begin(phi.matrix());
		const int iw = wp_;

		auto norm2 = [&](int ist) {
			auto res = gpu::run(1, gpu::reduce(n_pts), 0.0,
				[dV, mat_, ist] GPU_LAMBDA (auto, auto ip) {
					auto c = mat_[ip][ist];
					return dV*(inq::real(c)*inq::real(c) + inq::imag(c)*inq::imag(c));
				});
			INQKIT_GPU_SYNC();
			return res[0];
		};
		const double n2_before = norm2(iw);

		// ---- kick: psi_wp *= exp(+i dt_eff v(r)) --------------------------------
		// exp(-i dt (H - V)) ~ exp(-i dt H) exp(+i dt V): the correction kick has
		// the OPPOSITE sign to a potential's usual evolution phase because it
		// REMOVES v_SIC from the Hamiltonian. A real v leaves |psi|^2 untouched.
		{
			auto v_ = v.linear().cbegin();
			const double dt = dt_eff;
			gpu::run(n_pts, [dt, v_, mat_, iw] GPU_LAMBDA (auto ip) {
				const double ph = dt * v_[ip];
				const double c = cos(ph), s = sin(ph);
				auto w = mat_[ip][iw];
				const double wr = inq::real(w), wi = inq::imag(w);
				mat_[ip][iw] = inq::complex(wr*c - wi*s, wr*s + wi*c);
			});
			INQKIT_GPU_SYNC();
		}
		rep.kicked = true;

		// ---- Q projection: subtract every occupied component --------------------
		// Same modified-Gram-Schmidt pattern as WavePacket's injector. The
		// overlaps measured here are the leak the kick created (O(dt*v) per
		// step); their max is the gate observable max_overlap_pre.
		{
			auto const & occs = el.occupations()[0];
			const int n_states = phi.set_part().local_size();
			for(int j = 0; j < n_states; ++j) {
				if(j == iw) continue;
				if(occs[j] <= occ_thr_) continue;
				++rep.n_projected;

				auto res_re = gpu::run(1, gpu::reduce(n_pts), 0.0,
					[dV, mat_, j, iw] GPU_LAMBDA (auto, auto ip) {
						auto vj = mat_[ip][j]; auto vw = mat_[ip][iw];
						return dV*(inq::real(vj)*inq::real(vw) + inq::imag(vj)*inq::imag(vw));
					});
				INQKIT_GPU_SYNC();
				auto res_im = gpu::run(1, gpu::reduce(n_pts), 0.0,
					[dV, mat_, j, iw] GPU_LAMBDA (auto, auto ip) {
						auto vj = mat_[ip][j]; auto vw = mat_[ip][iw];
						return dV*(inq::real(vj)*inq::imag(vw) - inq::imag(vj)*inq::real(vw));
					});
				INQKIT_GPU_SYNC();
				const double re = res_re[0], im = res_im[0];
				rep.max_overlap_pre = std::max(rep.max_overlap_pre, std::sqrt(re*re + im*im));

				gpu::run(n_pts, [mat_, j, iw, re, im] GPU_LAMBDA (auto ip) {
					auto vj = mat_[ip][j]; auto vw = mat_[ip][iw];
					mat_[ip][iw] = inq::complex(
						inq::real(vw) - (re*inq::real(vj) - im*inq::imag(vj)),
						inq::imag(vw) - (re*inq::imag(vj) + im*inq::real(vj)));
				});
				INQKIT_GPU_SYNC();
			}
		}

		// ---- renormalise back to the pre-kick norm ------------------------------
		// The removed weight is REPORTED, never silently absorbed: if it stops
		// being << 1 the scheme is invalid and the gate catches it.
		{
			const double n2_after = norm2(iw);
			rep.norm_removed = (n2_before > 0.0) ? (1.0 - n2_after/n2_before) : 0.0;
			const double scale = std::sqrt(n2_before / n2_after);
			gpu::run(n_pts, [mat_, iw, scale] GPU_LAMBDA (auto ip) {
				auto w = mat_[ip][iw];
				mat_[ip][iw] = inq::complex(scale*inq::real(w), scale*inq::imag(w));
			});
			INQKIT_GPU_SYNC();
		}
		return rep;
	}

	// E_corrected assembly, variant-aware (plan §2 energy bookkeeping).
	double corrected_energy(double e_ks_total, SICStepReport const & r) const {
		switch(mode_) {
		case Mode::none:    return e_ks_total;
		case Mode::hartree: return e_ks_total - r.u_self;
		default:            return e_ks_total - r.u_self - r.exc_self;
		}
	}

	// -- internal helpers below. PUBLIC because nvcc forbids extended __device__
	// lambdas inside private member functions; the trailing underscore marks
	// them as implementation detail.

	// occ-weighted single-orbital density (mirrors jellium::orbital_density_field;
	// duplicated so the wavepacket module stays self-contained).
	template <class Electrons>
	inq::basis::field<inq::basis::real_space, double>
	wp_density_(Electrons const & el) const {
		auto const & phi = el.kpin()[0];
		inq::basis::field<inq::basis::real_space, double> n(phi.basis());
		n.fill(0.0);
		gpu::run(phi.basis().part().local_size(),
			[idx = wp_, occ = el.occupations()[0].cbegin(),
			 ph = phi.matrix().cbegin(), nn = n.linear().begin()] GPU_LAMBDA (auto ip) {
				nn[ip] = occ[idx] * norm(ph[ip][idx]);
			});
		return n;
	}

	// E_xc^unpol[n] alone (no potential needed).
	double exc_only_(inq::basis::field<inq::basis::real_space, double> const & n) const {
		inq::basis::field_set<inq::basis::real_space, double> dummy(n.basis(), 1);
		double exc = 0.0;
		eval_lda_unpol_(n, exc, dummy, /*want_v=*/false);
		return exc;
	}

	// v += v_xc^unpol[n]
	void add_vxc_unpol_(inq::basis::field<inq::basis::real_space, double> const & n,
	                    inq::basis::field<inq::basis::real_space, double> & v) const {
		inq::basis::field_set<inq::basis::real_space, double> vxc(n.basis(), 1);
		double exc = 0.0;
		eval_lda_unpol_(n, exc, vxc, /*want_v=*/true);
		gpu::run(n.basis().local_size(),
			[v_ = v.linear().begin(), x_ = vxc.matrix().cbegin()] GPU_LAMBDA (auto ip) {
				v_[ip] += x_[ip][0];
			});
		INQKIT_GPU_SYNC();
	}

	// Unpolarised LDA_X + LDA_C_PZ on a bare density, through INQ's OWN libxc
	// wrapper (hamiltonian::xc_term::evaluate_functional), so the subtraction
	// matches the run's Hamiltonian bit-for-bit (plan §0/D1). LDA needs no
	// gradient/laplacian/tau; the empty optionals mirror xc_term::operator().
	void eval_lda_unpol_(inq::basis::field<inq::basis::real_space, double> const & n,
	                     double & exc_out,
	                     inq::basis::field_set<inq::basis::real_space, double> & vxc_out,
	                     bool want_v) const {
		namespace ham = inq::hamiltonian;
		inq::basis::field_set<inq::basis::real_space, double> nset(n.basis(), 1);
		gpu::run(n.basis().local_size(),
			[m = nset.matrix().begin(), s = n.linear().cbegin()] GPU_LAMBDA (auto ip) {
				m[ip][0] = (s[ip] > 0.0) ? s[ip] : 0.0;
			});
		INQKIT_GPU_SYNC();

		auto grad = std::optional<decltype(inq::operations::gradient(nset))>{};
		auto lapl = std::optional<decltype(inq::operations::laplacian(nset))>{};
		auto ked  = std::optional<inq::basis::field<inq::basis::real_space, double>>{};
		auto vtau = std::optional<inq::basis::field_set<inq::basis::real_space, double>>{};

		inq::basis::field_set<inq::basis::real_space, double> vfunc(nset.skeleton());
		exc_out = 0.0;
		if(want_v) vxc_out.fill(0.0);

		for(int fid : {XC_LDA_X, XC_LDA_C_PZ}) {
			ham::xc_functional func(fid, 1);
			double e = 0.0;
			ham::xc_term::evaluate_functional(func, nset, grad, lapl, ked, e, vfunc, vtau);
			exc_out += e;
			if(want_v) {
				gpu::run(n.basis().local_size(),
					[o = vxc_out.matrix().begin(), f = vfunc.matrix().cbegin()] GPU_LAMBDA (auto ip) {
						o[ip][0] += f[ip][0];
					});
				INQKIT_GPU_SYNC();
			}
		}
	}

private:
	Mode mode_;
	int wp_;
	double occ_thr_;
};

} // namespace inqkit
