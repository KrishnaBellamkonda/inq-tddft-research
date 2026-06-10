/*
 * This file enumerates the degenerate orbital shells of a free-electron
 * jellium ground state in a cubic periodic box of side L, and selects a
 * fixed subset of "proxy" orbitals per shell for use in overlap and response
 * calculations.
 *
 * Shell structure
 * ---------------
 * Free electrons in an L³ box with periodic boundary conditions are
 * eigenstates of the kinetic operator indexed by integer triples (nₓ, nᵧ, n_z)
 * with eigenvalue:
 *
 *   E(n) = (2π/L)² (nₓ² + nᵧ² + n_z²) / 2
 *
 * Orbitals sharing the same |G|² = nₓ² + nᵧ² + n_z² are exactly degenerate
 * and form a shell. INQ stores Kohn-Sham orbitals sorted by eigenvalue, so
 * for a converged jellium ground state the INQ orbital index maps
 * deterministically to a shell: indices [0, 1) form shell 0 (|G|²=0),
 * indices [1, 7) form shell 1 (|G|²=1), and so on. The cumulative electron
 * counts at shell closures give the jellium magic numbers:
 *
 *   |G|²   shell size   cumulative electrons
 *      0            1              2
 *      1            6             14
 *      2           12             38
 *      3            8             54
 *      4            6             66
 *      5           24            114
 *      6           24            162   ← closed shell
 *
 * Not all integers appear as valid |G|² values. By Legendre's three-square
 * theorem, values of the form 4^k(8m+7) — such as 7, 15, 23, 28, 31, ... —
 * cannot be expressed as a sum of three integer squares and are absent from
 * the shell table.
 *
 * Proxy selection
 * ---------------
 * Within each shell, the first two orbitals (by INQ index) are selected as
 * proxies. The rationale is that in a jellium ground state the bath response
 * perturbation is approximately shell-symmetric, so any two members of a
 * degenerate shell carry similar information about the shell's collective
 * response. Picking two rather than one provides a within-shell consistency
 * check at modest computational cost. Proxies are always the lowest-index
 * members of their shell, making the selection deterministic and reproducible
 * across runs.
 *
 * Partial shells
 * --------------
 * When n_states does not coincide with a shell closure, the last shell is
 * truncated: its degeneracy and member list are clipped to the available
 * orbitals. This allows the same enumeration logic to work for any state
 * count, not only jellium magic numbers.
 *
 * CSV output
 * ----------
 * write_shells_csv() writes shells.csv alongside the overlap_proxies/
 * snapshots directory so that post-processing scripts can recover the
 * mapping from proxy column index back to (shell_id, gsq, degeneracy):
 *
 *   shell_id, gsq, degeneracy, n_proxies, proxy_indices,
 *   member_indices_first, member_indices_last
 *
 * API
 * ---
 *   // Enumerate all shells that fit within 101 orbitals:
 *   auto shells = inqkit::jellium::shells::enumerate_for_n_states(101);
 *
 *   // Select proxy orbitals (2 per shell by default):
 *   auto proxies = inqkit::jellium::shells::pick_proxies(shells);
 *
 *   // Write the shell/proxy mapping for post-processing:
 *   inqkit::jellium::shells::write_shells_csv(shells, proxies, output_dir);
 */
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
