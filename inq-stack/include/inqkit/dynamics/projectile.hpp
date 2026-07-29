/*
 * inqkit::dynamics::Projectile — a classical point/Gaussian projectile advanced by
 * Ehrenfest dynamics (velocity-Verlet) from the force it feels.
 *
 * WHY (twin-run Rung 2, docs/plans/twin-run-rung2-dynamic-spec.md): the classical
 * twin of a moving wavepacket. We do NOT drive it along a prescribed path — it
 * moves under its OWN Hellmann-Feynman force, so both twins evolve from identical
 * initial conditions and their divergence IS the quantum effect. The projectile is
 * a rigid Gaussian CHARGE realised as a moving `gaussian_projectile_perturbation`
 * (no ghost UPF, no r_cut aliasing); the force is computed wrapper-side each step
 * (see projectile_force.hpp) and fed here.
 *
 * This header is PURE (only inqkit::detail::Vec3, <cmath>) so the integrator is
 * unit-tested in isolation on the host — no INQ, no GPU. The force evaluation
 * (which needs INQ fields) lives in projectile_force.hpp.
 *
 * Integrator: velocity-Verlet in kick-drift-kick form, velocities synced at
 * integer steps (so KE = 1/2 m V^2 is the step-n kinetic energy). Symplectic:
 * energy is bounded, no secular drift. For a rigid symmetric Gaussian the self-
 * force vanishes by symmetry, so no self-force subtraction is needed.
 *
 * Callback contract (per RT step n, after the step, density reflects center R_n):
 *     Vec3 pos_n = proj.R();              // record proj_z = R_n (center used this step)
 *     Vec3 F     = projectile_force(...); // HF force at R_n on density_n
 *     proj.advance(F, dt);                // V -> V_n (completion kick), R -> R_{n+1}
 *     // record proj_vz = proj.V().z, ke = proj.ke()  (both at time n*dt)
 *     // moving perturbation now reads proj.R() = R_{n+1} for step n+1
 *
 * No inq/ or inq-study/ edit — wrapper-only.
 */
#ifndef INQKIT__DYNAMICS__PROJECTILE
#define INQKIT__DYNAMICS__PROJECTILE

#include <inqkit/detail/vec3.hpp>

namespace inqkit {
namespace dynamics {

using inqkit::detail::Vec3;

class Projectile {
public:
	// mass, charge in atomic units; R0, V0 initial position/velocity (Bohr, a.u.).
	Projectile(double mass, double charge, Vec3 R0, Vec3 V0)
		: mass_(mass), charge_(charge), R_(R0), V_(V0) {}

	// Advance one step given the force at the CURRENT position on the current
	// density (velocity-Verlet, KDK). After the call: V is the integer-step
	// velocity V_n, R is the next-step position R_{n+1}.
	void advance(Vec3 force, double dt) {
		Vec3 a = (1.0 / mass_) * force;
		if(started_) V_ += 0.5 * (a_last_ + a) * dt;   // complete previous drift's kick
		else         started_ = true;                   // first call: seed only
		a_last_ = a;
		R_ += V_ * dt + (0.5 * dt * dt) * a;            // drift to R_{n+1}
	}

	Vec3   R() const { return R_; }
	Vec3   V() const { return V_; }
	double mass() const { return mass_; }
	double charge() const { return charge_; }
	double ke() const { return 0.5 * mass_ * V_.norm2(); }

private:
	double mass_;
	double charge_;
	Vec3   R_{};
	Vec3   V_{};
	Vec3   a_last_{};
	bool   started_ = false;
};

} // namespace dynamics
} // namespace inqkit

#endif
