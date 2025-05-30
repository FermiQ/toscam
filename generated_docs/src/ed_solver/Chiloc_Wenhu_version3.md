# `src/ed_solver/Chiloc_Wenhu_version3.cpp`

## Overview

*   **Purpose:** Implements a revised version (v3) of the 'Chiloc' functionality, designed for computing local four-operator correlation functions (susceptibilities). This version features a significantly modified C-Fortran interface compared to `Chiloc_Wenhu_version1.cpp` and `Chiloc_Wenhu_version2.cpp`, particularly in how scalar arguments (passed by C++ reference) and the dimensions of matrix element arrays (passed explicitly) are handled. It also introduces a ground state energy parameter for relative energy calculations in the Boltzmann factor.
*   **Role in Project:** Provides a C++ computational kernel for an Exact Diagonalization (ED) solver. It takes many-body eigenvalues, matrix elements of creation operators, and now explicit total sizes for these matrix element arrays as input. It then computes specified components of the local susceptibility (\(\chi_{loc}\)) at given Matsubara frequencies.

## Key Components

*   **`inline static complex<double> PhiM_ii(double beta, double Ei, double Ej, double Ek, double El, complex<double> w1, complex<double> w2, complex<double> w3, double PHI_EPS)`**:
    Identical in implementation to versions 1 and 2. Calculates a specific term of a four-operator correlation function.
*   **`inline static complex<double> PhiM_ji(...)`**: Identical to versions 1 and 2.
*   **`inline static complex<double> PhiM_ki(...)`**: Identical to versions 1 and 2.
*   **`inline static complex<double> PhiM_li(...)`**: Identical to versions 1 and 2.
    *   *Overall purpose of `PhiM_*` functions*: These helper functions provide parts of the analytical formula for a three-frequency, four-operator correlation function, essential for building the full susceptibility.

*   **`extern "C" void chi_tilde_loc_(int &op, complex<double> &w1, complex<double> &w2, complex<double> &w3, double &PHI_EPS, double &beta, double &Z, double &gsE, int &sites, int &nup, int &ndn, double *cp_i_E, int &dim_E_i, ..., int &n1, ..., int &n28, double *cp_i_cdup, ..., complex<double> &pDNDN, complex<double> &pUPDN)`**:
    The main C++ function, exposed with C linkage for Fortran interoperability (name `chi_tilde_loc_`).
    *   **Key Interface Differences from v1 & v2:**
        *   **Scalar Arguments by Reference:** Most scalar input parameters (e.g., `op`, `w1`, `PHI_EPS`, `beta`, `Z`, `gsE`, `sites`, `nup`, `ndn`, and all `dim_E_*` which are dimensions of energy arrays) are now passed by C++ reference (`&`).
        *   **Ground State Energy `gsE`:** A new `double &gsE` parameter is introduced. The Boltzmann factor calculation is modified to `exp(-(cp_i_E[stati]-gsE)*beta)/Z`, making energies relative to `gsE`.
        *   **Array Pointers `double *`:** Energy and matrix element arrays are passed as flat `double *` pointers (similar to version 1, different from `double **` in version 2).
        *   **Explicit Matrix Array Sizes `n1...n28`:** A series of new integer arguments (`int &n1` through `int &n28`) are passed by reference. These explicitly provide the total number of elements for each corresponding 1D matrix element array (e.g., the array `cp_i_cdup` is assumed to have `n1*n2` elements, `cp_i_cddn` has `n3*n4` elements, and so on, as indicated by the commented-out debugging code). This is a major change, as previously these total sizes were calculated implicitly within the C++ code based on `dim_E_*` values.
    *   **Inputs (selected):**
        *   `op`: Integer flag (by ref) for selecting DNDN and/or UPDN calculation.
        *   `w1, w2, w3`: Complex Matsubara frequencies (by ref).
        *   `PHI_EPS, beta, Z, gsE`: Physical parameters (by ref).
        *   `sites, nup, ndn`: System/sector specification (by ref).
        *   `cp_i_E, dim_E_i, ...`: Pointers to energy arrays (`double*`) and their respective dimensions (`int&`).
        *   `n1, ..., n28`: Explicit total sizes for the flattened matrix element arrays (by ref).
        *   `cp_i_cdup, ...`: Pointers to flattened matrix element arrays (`double*`).
    *   **Outputs:**
        *   `pDNDN, pUPDN`: Calculated susceptibility components (by ref).
    *   **Functionality:** The core summation logic using `PhiM_*` functions and products of matrix elements appears structurally identical to previous versions. The main changes are in the function signature and how parameters (especially array dimensions) are received and how energies are treated in the Boltzmann factor.
    *   **Debugging Code:** The source contains extensive commented-out `std::cout` and `std::ofstream` statements, which were used for printing and verifying the passed arguments and intermediate values.

## Important Variables/Constants

*   **`double PHI_EPS` (passed by reference):** A small constant used in `PhiM_*` functions to avoid division by zero when frequency/energy denominators are close to zero.
*   **`double gsE` (passed by reference):** The ground state energy of the system. Used to calculate energies relative to the ground state for the Boltzmann factor `exp(-(cp_i_E[stati]-gsE)*beta)/Z`.
*   **`int n1, ..., n28` (passed by reference):** These are new crucial inputs in this version. They represent the total pre-calculated sizes of the 1D arrays that store the matrix elements. For example, if `cp_i_cdup` stores a matrix of dimensions `dim_destination_sector` x `dim_source_sector`, then `n1*n2` (or whichever variables correspond to `cp_i_cdup`) would be `dim_destination_sector * dim_source_sector`. The Fortran caller is responsible for providing these correct total sizes.
*   Scalar arguments being passed by C++ reference is a key characteristic of this version's interface.
*   Energy and matrix element data arrays are passed as flat `double *` pointers.

## Usage Examples

This C++ function `chi_tilde_loc_` is designed to be called from a Fortran ED solver.

**Key considerations for Fortran integration:**

1.  The Fortran ED solver calculates many-body eigenstates, eigenvalues, and matrix elements. It also needs to determine the ground state energy `gsE`.
2.  Scalar arguments (integers, doubles, complex values) should be passed directly from Fortran; the C++ interface will handle them as references (effectively, Fortran passes their memory addresses).
3.  Energy arrays and matrix element arrays should be passed as standard Fortran arrays (which are contiguous and will appear as `double *` to C++).
4.  Crucially, Fortran must pre-calculate and pass the total number of elements for each matrix element array through the `n1, ..., n28` arguments. For example, if `cp_i_cdup_fortran` is a 2D Fortran array `DIM1 x DIM2`, it would be passed as `cp_i_cdup_fortran(1,1)` (the starting address), and the corresponding `n_row_arg`, `n_col_arg` (representing `n1`, `n2` for this specific matrix) would be `DIM1`, `DIM2`. The C++ side then uses `n_row_arg * n_col_arg` as the total size for its 1D access. *(Correction based on C++ code: The C++ code still uses 1D indexing like `cp_i_cddn[dimi*j+stati]`. The `n1..n28` parameters are more likely the explicit total sizes of these flattened arrays, not individual dimensions to be multiplied in C++.)* The commented-out C++ debug lines like `file<<"cp_i_cdup="<<endl; for (int i=0; i<n1*n2; i++ )` confirm that `n1*n2` is the total size of the `cp_i_cdup` 1D array.

**Conceptual C++ side of the `extern "C"` function:**
```cpp
// extern "C" void chi_tilde_loc_(
//     int &op, std::complex<double> &w1, /* ... other scalar args by ref ... */ double &gsE,
//     double *cp_i_E, int &dim_E_i, /* ... */
//     int &n1_mult_n2, /* effectively the total size for cp_i_cdup */
//     /* ... up to n27_mult_n28 for cp_m2dn_cddn ... */
//     double *cp_i_cdup, /* ... other matrix element arrays as double* ... */
//     std::complex<double> &pDNDN, std::complex<double> &pUPDN) {
//
//     // Access scalar as: op, beta, Z, gsE (these are references to Fortran variables)
//     // Access energy array as: cp_i_E[stati]
//     // Access matrix element array cp_i_cdup (which has total size n1_mult_n2 from Fortran)
//     // as: cp_i_cdup[some_1D_index_up_to_n1_mult_n2-1]
// }
```
Developers need to ensure the Fortran call correctly matches this modified interface, especially the passing of scalars by reference and providing the correct total sizes for the flattened matrix element arrays.

## Dependencies and Interactions

*   **Internal Dependencies:**
    *   The `PhiM_*` helper functions are defined and used within this same file.

*   **External Libraries:**
    *   **`<complex>`**: Standard C++ header for `std::complex<double>`.
    *   **`<iostream>`, `<fstream>`, `<iomanip>`**: Included for I/O operations, primarily used in the (now commented-out) debugging sections of the code for printing values to console or files.
    *   **`using namespace std;`**: Standard practice for using names from the C++ standard library.

*   **Interactions with other components:**
    *   **Fortran ED Solver (Primary Caller):** `chi_tilde_loc_` is intended to be called from Fortran.
    *   **Data Flow:**
        1.  The Fortran solver computes eigenstates, eigenvalues (including `gsE`), and matrix elements.
        2.  It passes scalar parameters by address (handled by C++ references). Energy and matrix element arrays are passed as flat `double*`. Crucially, the total sizes of these flattened matrix element arrays are passed via `n1, ..., n28`.
        3.  The C++ function calculates `pDNDN` and `pUPDN`.
        4.  Results are returned to Fortran through the reference parameters `pDNDN` and `pUPDN`.
    *   **Key Interface Changes from v1/v2:**
        *   Scalar arguments are passed by C++ reference.
        *   `gsE` is a new input, modifying the Boltzmann factor.
        *   Matrix element arrays are `double*` (like v1, unlike v2's `double**`).
        *   Explicit total sizes (`n1..n28`) for the 1D representation of matrix element arrays are now required inputs, replacing implicit calculation of these sizes based on `dim_E_*` values.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
