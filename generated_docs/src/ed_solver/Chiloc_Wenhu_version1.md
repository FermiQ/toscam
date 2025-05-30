# `src/ed_solver/Chiloc_Wenhu_version1.cpp`

## Overview

*   **Purpose:** Likely implements a specific version or part of the 'Chiloc' functionality, possibly related to vertex calculations or a particular solution strategy for computing local correlation functions (like a two-particle Green's function or susceptibility) in an exact diagonalization context. The name 'Wenhu_version1' suggests it's a contribution or specific variant by Wenhu. The detailed comments suggest it calculates \(\Phi(E_i, E_j, E_k, E_l; \omega_1, \omega_2, \omega_3)\), which is a Fourier transform of a four-operator expectation value \(\langle T c c c c \rangle\).
*   **Role in Project:** Contributes a specialized C++ component to the Exact Diagonalization solver suite, focused on calculating specific four-operator correlation functions (Chi_loc) based on eigenstates and eigenvalues obtained from an ED solver.

## Key Components

*   **`inline static complex<double> PhiM_ii(double beta, double Ei, double Ej, double Ek, double El, complex<double> w1, complex<double> w2, complex<double> w3, double PHI_EPS)`**:
    Calculates a specific term of a four-operator correlation function for Matsubara frequencies. It takes the inverse temperature `beta`, energies of four many-body states (`Ei, Ej, Ek, El`), three complex Matsubara frequencies (`w1, w2, w3`), and a small epsilon `PHI_EPS` to handle potential divisions by zero. The `_ii` suffix likely denotes a specific combination or ordering of operators/indices in the full expression.
*   **`inline static complex<double> PhiM_ji(...)`**: Similar to `PhiM_ii`, calculates a different term (denoted by `_ji`) for the same set of input parameters.
*   **`inline static complex<double> PhiM_ki(...)`**: Similar to `PhiM_ii`, calculates a different term (denoted by `_ki`).
*   **`inline static complex<double> PhiM_li(...)`**: Similar to `PhiM_ii`, calculates a different term (denoted by `_li`).
    *   *Overall purpose of `PhiM_*` functions*: These functions implement the analytical formula for different parts of a Lehmann representation of a three-frequency, four-operator correlation function, handling cases where energy denominators might become small.

*   **`void chi_tilde_loc(int op, complex<double> w1, complex<double> w2, complex<double> w3, double PHI_EPS, double beta, double Z, int sites, int nup, int ndn, double *cp_i_E, int dim_E_i, double *cp_pup_E, int dim_E_pup, ..., complex<double> &pDNDN, complex<double> &pUPDN)`**:
    This is the main function in the file. It calculates components of the local susceptibility, specifically \(\chi_{DNDN}\) (density-density) and/or \(\chi_{UPDN}\) (spin-spin like), based on the input `op` flag.
    *   **Inputs:**
        *   `op`: Integer flag to select which correlation functions to compute (bitwise: `&1` for DNDN, `&2` for UPDN).
        *   `w1, w2, w3`: Three complex Matsubara frequencies.
        *   `PHI_EPS`: Small constant to avoid division by zero in `PhiM_*` calls.
        *   `beta`: Inverse temperature.
        *   `Z`: Partition function.
        *   `sites`: Number of sites in the system (likely for context, not directly used in core math).
        *   `nup, ndn`: Number of up and down electrons defining the initial sector for state `i`.
        *   `cp_i_E, dim_E_i, ...`: A series of pointers to arrays containing energy eigenvalues for different electronic sectors (e.g., `cp_i_E` for sector `(nup, ndn)`, `cp_pup_E` for `(nup+1, ndn)`, `cp_p2dn_E` for `(nup, ndn+2)`, etc.), along with the dimensions of these energy lists.
        *   `cp_i_cdup, cp_i_cddn, ...`: A series of pointers to arrays containing matrix elements of creation operators (e.g., `cp_i_cdup` for \(\langle \text{nup+1,ndn} | c^\dagger_{up} | \text{nup,ndn} \rangle\)). These are stored as 1D arrays, and the indexing convention is described in comments (e.g., `cp_i_cdup[dimj*i+j]`).
    *   **Outputs:**
        *   `pDNDN`: Result for the density-density like susceptibility, passed by reference.
        *   `pUPDN`: Result for the up-down (spin-spin like) susceptibility, passed by reference.
    *   **Functionality:**
        The function iterates through all initial many-body states `stati` in the sector `(nup, ndn)`. For each `stati`, it calculates contributions to `cDNDN` and `cUPDN`. This involves nested loops over intermediate states `j, k, l` belonging to various other sectors reachable by applying creation/destruction operators. Inside these loops, it calls the appropriate `PhiM_*` functions with the corresponding energies and frequencies and multiplies their results by products of the pre-calculated matrix elements of creation operators. The structure of these summations, indicated by comments like `// tot =0/20`, `// tot =2/22`, corresponds to summing over different terms in the full expression of the susceptibility, likely arising from different Wick contractions or time orderings of the four operators. Finally, it performs a Boltzmann weighting (`exp(-cp_i_E[stati]*beta)/Z`) for each initial state's contribution.

## Important Variables/Constants

*   **`double PHI_EPS`**: This is a crucial input parameter for `chi_tilde_loc` and subsequently for all `PhiM_*` functions. It's a small positive double value used as a threshold to determine if an energy denominator `w + E_state_diff` is effectively zero. If `abs(w + E_state_diff) <= PHI_EPS`, alternative formulas (often involving `beta` or squared terms) are used to handle these degenerate or near-degenerate cases, preventing division by zero and correctly capturing terms that arise from such degeneracies in perturbation theory.
*   Other variables are local to the functions or passed as arguments. The core data consists of the many-body energy eigenvalues and matrix elements of creation operators, which are provided by the calling (Fortran) code.

## Usage Examples

This C++ code is intended to be compiled and linked with a Fortran-based Exact Diagonalization (ED) solver. The Fortran code would call the `chi_tilde_loc` function, likely through a C-Fortran interface layer.

**Conceptual Fortran call sequence:**

1.  The Fortran ED solver calculates many-body eigenstates and eigenvalues for all relevant Fock sectors (e.g., `(Nup,Ndn)`, `(Nup+1,Ndn)`, `(Nup,Ndn+1)`, etc.).
2.  It then computes matrix elements of creation operators (e.g., \( \langle \Psi_{sectorA, stateA} | c^\dagger_{\sigma, site_m} | \Psi_{sectorB, stateB} \rangle \)) between these eigenstates. These are stored in arrays.
3.  The Fortran code sets up the required parameters: `op` flag, Matsubara frequencies `w1, w2, w3`, `PHI_EPS`, `beta`, `Z`, `sites`, `nup`, `ndn` for the initial sector.
4.  It then calls a wrapper function (interfacing Fortran and C++) that ultimately invokes `chi_tilde_loc`, passing pointers to the energy arrays and matrix element arrays.

**Example C++ wrapper for Fortran (illustrative, not in the provided file):**
```cpp
// extern "C" void compute_chiloc_wenhu_version1_wrapper(
//     int op, double w1r, double w1i, double w2r, double w2i, double w3r, double w3i,
//     double PHI_EPS, double beta, double Z, int sites, int nup, int ndn,
//     double *cp_i_E, int dim_E_i, /* ... many more energy and matrix element arrays ... */
//     double* out_DNDN_real, double* out_DNDN_imag,
//     double* out_UPDN_real, double* out_UPDN_imag) {
//
//     std::complex<double> w1_c(w1r, w1i);
//     std::complex<double> w2_c(w2r, w2i);
//     std::complex<double> w3_c(w3r, w3i);
//     std::complex<double> pDNDN_c, pUPDN_c;
//
//     // Note: The order and number of cp_*_E and cp_*_c* arrays must match the chi_tilde_loc signature exactly.
//     chi_tilde_loc(op, w1_c, w2_c, w3_c, PHI_EPS, beta, Z, sites, nup, ndn,
//                   cp_i_E, dim_E_i, /* ... all other required array pointers and dimensions ... */
//                   pDNDN_c, pUPDN_c);
//
//     *out_DNDN_real = pDNDN_c.real(); *out_DNDN_imag = pDNDN_c.imag();
//     *out_UPDN_real = pUPDN_c.real(); *out_UPDN_imag = pUPDN_c.imag();
// }
```
The Fortran code would then receive the real and imaginary parts of `pDNDN` and `pUPDN`.
Further examples specific to its integration and the exact structure of the Fortran call are to be added by developers.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   The file is self-contained in terms of C++ code (the `PhiM_*` functions are used by `chi_tilde_loc` within the same file). No custom project-specific headers (`#include "my_header.h"`) are present.

*   **External Libraries:**
    *   **`<complex>`**: This standard C++ header is included for using `std::complex<double>` to represent complex numbers (frequencies, matrix elements, and results).
    *   **`using namespace std;`**: This directive brings names from the standard namespace into the global scope, implying potential use of other standard library components like `abs()`, `exp()` (though these are often available for `double` via `<cmath>` which might be implicitly included or not needed if only `std::complex` methods are used).

*   **Interactions with other components:**
    *   **Fortran ED Solver (Primary Caller):** This C++ code is designed to be a computational kernel called from a larger, likely Fortran-based, Exact Diagonalization (ED) program.
    *   **Data Flow:**
        1.  The Fortran ED solver performs the initial diagonalization, obtaining many-body eigenstates and their energies for various Fock space sectors.
        2.  It then computes matrix elements of fermionic creation (and implicitly, destruction) operators between these eigenstates.
        3.  These extensive datasets (energies and matrix elements, which can be large) are passed as pointers to arrays into the C++ function `chi_tilde_loc`.
        4.  Physical parameters like inverse temperature (`beta`), partition function (`Z`), and the Matsubara frequencies (`w1, w2, w3`) at which the susceptibility is to be computed are also passed from Fortran.
        5.  The C++ code computes the specified components of the local susceptibility (`pDNDN`, `pUPDN`).
        6.  These complex-valued results are passed back to the Fortran caller, typically by filling values at memory locations pointed to by output arguments.
    *   **Purpose of C++ Implementation:** This specific part of the calculation (summing many terms involving complex numbers and products of matrix elements) might be implemented in C++ for reasons such as:
        *   Easier handling of `std::complex` arithmetic compared to Fortran's complex types in some contexts.
        *   Potentially leveraging C++ features for performance or code organization for this specific, complex summation.
        *   Historical reasons or developer preference for this computational kernel.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
