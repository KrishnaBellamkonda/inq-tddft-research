// ============================================================================
// inqkit/jellium/shells.hpp
//
// Enumerate the degenerate orbital shells of a free-electron jellium ground
// state on a cubic box and pick a fixed set of "proxy" orbitals per shell.
//
// For free electrons in an L^3 cubic box with periodic BC, the eigenstates
// are plane waves indexed by integer triples (n_x, n_y, n_z) with eigenvalue
//
//     E(n)  =  (2 pi / L)^2  *  (n_x^2 + n_y^2 + n_z^2)  /  2
//
// Degenerate shells are levels of constant |G|^2 = n_x^2 + n_y^2 + n_z^2.
// The cumulative shell structure for jellium magic numbers:
//
//     |G|^2  shell_size  cumulative spatial states  cumulative electrons
//        0             1                          1                     2
//        1             6                          7                    14
//        2            12                         19                    38
//        3             8                         27                    54
//        4             6                         33                    66
//        5            24                         57                   114
//        6            24                         81                   162    <-- N=162 closed shell
//        8             6                         87                   174
//        9            30                        117                   234
//
// The "magic numbers" 2, 14, 38, 54, 66, 114, 162 are closed shells.
// |G|^2 = 7 is impossible (no integer (n_x,n_y,n_z) satisfies it: Lagrange's
// 4-square theorem). Same for |G|^2 = 15, 23, 28, 31, ... (excluded levels).
//
// SHELL ASSIGNMENT IN INQ ORBITAL ORDER. INQ stores Kohn-Sham orbitals
// sorted by eigenvalue at the converged GS. For a uniform jellium GS this
// matches the |G|^2 ordering exactly (lowest |G|^2 first). The mapping from
// INQ orbital_index to shell_id is therefore deterministic given (L, N, n_states):
//   indices [0, 1)        -> shell 0
//   indices [1, 7)        -> shell 1   (size 6)
//   indices [7, 19)       -> shell 2   (size 12)
//   indices [19, 27)      -> shell 3   (size 8)
//   indices [27, 33)      -> shell 4   (size 6)
//   indices [33, 57)      -> shell 5   (size 24)
//   indices [57, 81)      -> shell 6   (size 24)
//   indices [81, 87)      -> shell 8   (size 6)   <-- shell 7 absent
//   indices [87, 117)     -> shell 9   (size 30)
//   ...
//
// Within a shell we pick the FIRST TWO orbitals as proxies. The user-supplied
// rationale: in jellium the bath response perturbation is approximately
// shell-symmetric, so any two members of a degenerate shell carry similar
// information about how the shell as a whole responds. Picking 2 (vs 1)
// gives a within-shell consistency check at modest cost.
//
// API:
//
//   auto shells = inqkit::jellium::shells::enumerate_for_n_states(101);
//   // shells is a vector<ShellInfo> with .id, .gsq, .degeneracy,
//   // .first_index, .members (full list of orbital indices in the shell).
//   // Truncated at the last shell that fits within n_states.
//
//   auto proxies = inqkit::jellium::shells::pick_proxies(shells);
//   // vector<int> of orbital indices, 2 per shell (or 1 for shell-of-1).
//
//   inqkit::jellium::shells::write_shells_csv(shells, proxies, dir);
//   // dir/shells.csv with columns:
//   //   shell_id, gsq, degeneracy, n_proxies, proxy_indices,
//   //   member_indices_first, member_indices_last
// ============================================================================
#pragma once

#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>
#include <filesystem>

namespace inqkit::jellium::shells {

struct ShellInfo {
    int shell_id;          // ordinal index 0, 1, 2, ... in eigenvalue order
    int gsq;               // |G|^2 = n_x^2 + n_y^2 + n_z^2
    int degeneracy;        // # of orbitals in the shell
    int first_index;       // first orbital_index in INQ's order
    std::vector<int> members;  // orbital_indices [first_index, first_index+degeneracy)
};

// Hard-coded shell table for jellium in a cubic box.
// Pairs of (|G|^2, degeneracy) up to a comfortably high level. The table
// covers everything up to the v2-run state count (101) plus headroom for
// future denser runs.
//
// Excluded |G|^2 (no integer triple): 7, 15, 23, 28, 31, 39, 47, 55, 60, ...
// (Legendre's three-square theorem: |G|^2 = 4^k * (8m+7) is excluded.)
//
inline std::vector<std::pair<int,int>> default_shell_table() {
    return {
        { 0,  1},
        { 1,  6},
        { 2, 12},
        { 3,  8},
        { 4,  6},
        { 5, 24},
        { 6, 24},
        { 8,  6},
        { 9, 30},
        {10, 24},
        {11, 24},
        {12,  8},
        {13, 24},
        {14, 48},
        {16,  6},
        {17, 48},
        {18, 36},
        {19, 24},
        {20, 24},
        {21, 48},
        {22, 24},
        {24, 24},
    };
}

// Enumerate shells until cumulative orbital count >= n_states. Truncate
// the last shell so the cumulative orbital count == n_states. (If
// truncation lands inside a shell, that shell appears in the output but
// its degeneracy is reduced and members are clipped.)
inline std::vector<ShellInfo>
enumerate_for_n_states(int n_states)
{
    std::vector<ShellInfo> out;
    auto table = default_shell_table();
    int cumulative = 0;
    for (std::size_t k = 0; k < table.size() && cumulative < n_states; ++k) {
        int gsq = table[k].first;
        int deg = table[k].second;
        int first = cumulative;
        int n_take = std::min(deg, n_states - cumulative);
        ShellInfo info;
        info.shell_id   = static_cast<int>(k);
        info.gsq        = gsq;
        info.degeneracy = n_take;
        info.first_index = first;
        info.members.resize(n_take);
        for (int m = 0; m < n_take; ++m)
            info.members[m] = first + m;
        out.push_back(std::move(info));
        cumulative += n_take;
    }
    if (cumulative < n_states) {
        throw std::runtime_error(
            "inqkit::jellium::shells: shell table exhausted before reaching "
            "n_states; extend default_shell_table()");
    }
    return out;
}

// Pick proxies: 2 per shell when degeneracy >= 2, else 1 (the only one).
// Proxies are the FIRST orbitals of each shell — fixed for reproducibility
// across runs.
inline std::vector<int>
pick_proxies(std::vector<ShellInfo> const& shells, int n_per_shell = 2)
{
    std::vector<int> out;
    for (auto const& s : shells) {
        int n_take = std::min(n_per_shell, s.degeneracy);
        for (int m = 0; m < n_take; ++m)
            out.push_back(s.first_index + m);
    }
    return out;
}

// Write shells.csv next to overlap_proxies/ snapshots so the postprocess
// can recover the (proxy column -> shell, degeneracy) mapping.
inline void
write_shells_csv(std::vector<ShellInfo> const& shells,
                 std::vector<int>       const& proxies,
                 std::string            const& dir)
{
    std::filesystem::create_directories(dir);
    std::ofstream f(dir + "/shells.csv");
    if (!f)
        throw std::runtime_error("write_shells_csv: cannot open "
                                 + dir + "/shells.csv");
    f << "shell_id,gsq,degeneracy,n_proxies,proxy_indices,"
         "member_indices_first,member_indices_last\n";

    // For each shell, find which proxies belong to it (members.first to
    // first+degeneracy), preserving the order in `proxies`.
    for (auto const& s : shells) {
        std::vector<int> sps;
        for (int p : proxies)
            if (p >= s.first_index && p < s.first_index + s.degeneracy)
                sps.push_back(p);
        f << s.shell_id << ","
          << s.gsq << ","
          << s.degeneracy << ","
          << sps.size() << ",\"";
        for (std::size_t i = 0; i < sps.size(); ++i) {
            f << sps[i];
            if (i + 1 < sps.size()) f << ",";
        }
        f << "\","
          << s.first_index << ","
          << (s.first_index + s.degeneracy - 1) << "\n";
    }
}

}  // namespace inqkit::jellium::shells
