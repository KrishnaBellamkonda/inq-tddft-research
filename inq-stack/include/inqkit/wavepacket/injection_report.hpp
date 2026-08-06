/*
 * A class that keeps track of the wave packet information together.
 * This object can be used to make sanity checks and ensure the wave
 * packet is behaving as expected. 
 * */

#pragma once

#include <cmath>
#include <vector>

namespace inqkit {

struct InjectionReport {
    int    kpoint_index     = 0;
    int    state_index      = -1;
    double norm_before      = 0.0;   // norm of the orbital slot before injection
    double norm_after       = 0.0;   // norm of WP after injection (and ortho if applied)
    double max_overlap      = 0.0;   // max |<psi_i|psi_wp>| over occupied states (pre-ortho)
    bool   orthogonalised   = false;
    bool   passed_tolerance = false;

    // ---- orthogonalisation LOSS -------------------------------------------
    // How much of the packet the Gram-Schmidt projection carved away. Needed
    // because norm_after is measured AFTER the renormalisation step and is
    // therefore ~1 by construction: it cannot express the deformation, and
    // max_overlap only reports the single LARGEST overlap, not the total.
    //
    // Launching a wavepacket close to a metal surface puts it inside the
    // electronic spill-out, where it acquires real overlap with the occupied
    // manifold; the projection then removes that component and the renormalise
    // hides the loss. removed_weight is the honest measure of "how much of the
    // Gaussian survived as a Gaussian".
    double norm_pre_ortho   = 0.0;   // ||psi|| of the raw Gaussian, before GS
    double norm_pre_renorm  = 0.0;   // ||psi|| after GS, BEFORE renormalisation
    double removed_weight   = 0.0;   // 1 - (norm_pre_renorm/norm_pre_ortho)^2
    double sum_overlap_sq   = 0.0;   // sum_i |<psi_i|psi_wp>|^2 on the FIRST pass

    // |<psi_i|psi_wp>| per state i, first pass, in state order. Says WHICH bath
    // states the packet is mixing with (e.g. surface vs bulk-like), which is the
    // diagnostic you want the moment removed_weight comes out too large.
    std::vector<double> overlap_by_state{};

    // Independent cross-check of removed_weight. The occupied KS states are
    // mutually orthonormal, so in exact arithmetic subtracting all projections
    // removes exactly sum_i |<psi_i|psi_wp>|^2 of the squared norm:
    //
    //     sum_overlap_sq  ==  norm_pre_ortho^2 - norm_pre_renorm^2
    //
    // The two sides are computed by completely different reductions, so their
    // agreement is a genuine known-case test that runs on every real injection.
    // Returns the ABSOLUTE discrepancy between the two routes.
    double ortho_closure_residual() const {
        return std::abs(sum_overlap_sq -
                        (norm_pre_ortho * norm_pre_ortho -
                         norm_pre_renorm * norm_pre_renorm));
    }
};

} // namespace inqkit
