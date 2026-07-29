/*
 * A class that keeps track of the wave packet information together.
 * This object can be used to make sanity checks and ensure the wave
 * packet is behaving as expected. 
 * */

#pragma once

namespace inqkit {

struct InjectionReport {
    int    kpoint_index     = 0;
    int    state_index      = -1;
    double norm_before      = 0.0;   // norm of the orbital slot before injection
    double norm_after       = 0.0;   // norm of WP after injection (and ortho if applied)
    double max_overlap      = 0.0;   // max |<psi_i|psi_wp>| over occupied states (pre-ortho)
    bool   orthogonalised   = false;
    bool   passed_tolerance = false;
};

} // namespace inqkit
