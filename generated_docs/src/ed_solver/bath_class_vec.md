# `src/ed_solver/bath_class_vec.f90`

## Overview

*   **Purpose:** This module provides utility subroutines to convert the parameters of a `bath_type` object (which are stored in various `masked_matrix_type` structures representing bath energies, hybridizations, and pairing terms) into a single, flat, 1D real vector. It also provides the reverse operation: updating the `bath_type` object's matrices from such a flat vector.
*   **Role in Project:** Serves as a crucial interface layer. Numerical optimization algorithms (e.g., conjugate gradient, L-BFGS, as might be used in fitting the bath hybridization function in `bath_class_hybrid`) typically operate on a simple 1D vector of real parameters. This module handles the complex task of "flattening" the structured and potentially complex-valued, hermitian matrices of the `bath_type` into such a vector, and then "unflattening" the vector back into the `bath_type` structure after optimization.

## Key Components

*   **`MODULE bath_class_vec`**:
    The module containing the public subroutines for parameter vectorization.

*   **Public Subroutines:**
    *   **`SUBROUTINE bath2vec(BATH)`**:
        Converts the parameters from the structured `masked_matrix_type` components within the input `BATH` object into its `BATH%vec` 1D `REAL(DP)` array.
        *   **Argument:** `BATH (TYPE(bath_type), INTENT(INOUT))`. The `BATH%vec` component is populated.
        *   **Functionality:**
            *   Initializes `BATH%vec` if not already allocated.
            *   Iterates through the relevant `masked_matrix_type` components of `BATH`:
                1.  `BATH%Vbc(spin=1)` (Impurity-Bath Hybridization, spin 1)
                2.  `BATH%PVbc(spin=1)` (Impurity-Bath Pairing Hybridization, spin 1)
                3.  `BATH%Eb(spin=1)` (Bath Energy, spin 1)
                4.  `BATH%Vbc(spin=2)`
                5.  `BATH%PVbc(spin=2)`
                6.  `BATH%Eb(spin=2)`
                7.  `BATH%Pb` (Bath Pairing, if `BATH%SUPER` is true)
            *   For each matrix, it calls the internal `add_mm2vec` subroutine to append its independent real parameters to `BATH%vec`.
            *   On the first call (controlled by `first_call` and `ALL_FIRST_CALL`), it performs a consistency check to ensure that the total number of parameters gathered (`offset`) matches the pre-calculated `BATH%nparam`.

    *   **`SUBROUTINE vec2bath(BATH)`**:
        Updates the parameters within the `masked_matrix_type` components of the input `BATH` object using values from its `BATH%vec` 1D real vector.
        *   **Argument:** `BATH (TYPE(bath_type), INTENT(INOUT))`. Its matrix components are updated from `BATH%vec`.
        *   **Functionality:**
            *   Iterates through the matrices in the same order as `bath2vec`.
            *   For each matrix, it calls the internal `extract_mmvec_from_vec` subroutine to read its parameters from the appropriate segment of `BATH%vec`.
            *   After extracting parameters for all defined matrices, it includes logic to fill missing spin components if one spin channel has parameters defined but the other does not (e.g., if `Vbc(2)` has no independent parameters but `Vbc(1)` does, it copies parameters from `Vbc(1)` to `Vbc(2)` using `vec2masked_matrix`). This implies an assumption of spin symmetry if one channel is unparameterized.
            *   Performs a consistency check on `BATH%nparam` on the first call.

*   **Internal Subroutines (Helper routines within the module):**
    *   **`SUBROUTINE add_mm2vec(offset, MM, vec)`**:
        Appends the independent parameters of a single `masked_matrix_type` (`MM`) to the 1D `vec`, updating `offset`.
        *   Calls `masked_matrix2vec(MM)` (from `masked_matrix_class`) to ensure `MM%rc%vec` (and `vecdiag`, `vecoffdiag` for complex hermitian) is current.
        *   Handles complex matrices (if `#ifdef _complex` is active):
            *   Non-Hermitian: Appends real parts of `MM%rc%vec`, then imaginary parts. `offset` increases by `2 * nind`.
            *   Hermitian: Appends real diagonal parts (`MM%rc%vecdiag`), then real parts of off-diagonal elements (`MM%rc%vecoffdiag`), then imaginary parts of off-diagonal elements. `offset` increases by `ninddiag + 2 * nindoffdiag`.
        *   Handles real matrices (if `_complex` is not defined): Appends all elements of `MM%rc%vec` directly. `offset` increases by `nind`.
    *   **`SUBROUTINE extract_mmvec_from_vec(offset, MM, vec)`**:
        Extracts parameters for a single `masked_matrix_type` (`MM`) from the 1D `vec`, starting at `offset`, and updates `MM`.
        *   Populates `MM%rc%vecdiag` and `MM%rc%vecoffdiag` (for complex hermitian) or `MM%rc%vec` (for other cases) from the appropriate segment of `vec`.
        *   Handles complex/real and hermitian/non-hermitian cases in a way that reverses the logic of `add_mm2vec`.
        *   Updates `offset`.
        *   Crucially, after populating the vector representations within `MM%rc`, it calls `vec2masked_matrix(MM)` (from `masked_matrix_class`) to reconstruct the full matrix `MM%rc%mat` from these vector components.
    *   **`SUBROUTINE write_bathvec(BATH, UNIT)`**: (Not public from the module scope as defined, but could be). Formats and writes the content of `BATH%vec` to a specified Fortran `UNIT` or to the default log unit.

## Important Variables/Constants

*   **`BATH%vec :: REAL(DP), POINTER`**: This is the central data array managed by this module. It's a component of `bath_type` that `bath2vec` fills from the structured matrix data, and `vec2bath` uses to repopulate the structured matrix data.
*   **`offset :: INTEGER`**: A local variable within `bath2vec`, `vec2bath`, `add_mm2vec`, and `extract_mmvec_from_vec`. It acts as a running index or pointer to the current position within the flat `BATH%vec` array, ensuring that parameters for different matrices are packed and unpacked correctly without overlap.
*   **`first_call :: LOGICAL, SAVE = .true.`**: A saved logical variable used in `bath2vec` and `vec2bath` to perform a consistency check between the calculated total number of parameters (`offset`) and `BATH%nparam` only during the very first invocation of these routines (or first after `ALL_FIRST_CALL` from `globalvar_ed_solver` resets it).

## Usage Examples

```fortran
MODULE example_bath_vec_usage
  USE bath_class, ONLY: bath_type, read_bath, delete_bath ! Assuming read_bath initializes nparam
  USE bath_class_vec
  USE genvar, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_bath_vectorization

CONTAINS
  SUBROUTINE test_bath_vectorization(bath_input_file, num_imp_sites)
    CHARACTER(LEN=*), INTENT(IN) :: bath_input_file
    INTEGER, INTENT(IN) :: num_imp_sites
    TYPE(bath_type) :: my_bath
    REAL(DP), ALLOCATABLE :: optimized_params_vec(:)
    INTEGER :: i

    ! 1. Initialize bath_type (e.g., from a file, which also sets my_bath%nparam)
    CALL read_bath(my_bath, bath_input_file, num_imp_sites)
    PRINT *, "Initial bath nparam: ", my_bath%nparam

    ! 2. Convert structured bath parameters into the flat vector my_bath%vec
    CALL bath2vec(my_bath)
    PRINT *, "Bath parameters flattened into my_bath%vec. Size: ", SIZE(my_bath%vec)
    IF (SIZE(my_bath%vec) > 0) THEN
      PRINT *, "First few elements of my_bath%vec: ", my_bath%vec(1:MIN(5, SIZE(my_bath%vec)))
    END IF

    ! 3. (Simulate optimization) A generic minimizer would operate on my_bath%vec.
    !    Let's assume the minimizer produces a new vector 'optimized_params_vec'.
    IF (ASSOCIATED(my_bath%vec)) THEN
      ALLOCATE(optimized_params_vec(SIZE(my_bath%vec)))
      ! For demonstration, just slightly modify the parameters
      optimized_params_vec = my_bath%vec * 1.1
      my_bath%vec = optimized_params_vec
      DEALLOCATE(optimized_params_vec)
    END IF

    ! 4. Update the structured bath_type (Eb, Vbc, Pb matrices) from the modified my_bath%vec
    CALL vec2bath(my_bath)
    PRINT *, "Bath matrices updated from modified my_bath%vec."
    ! Now my_bath's matrix components (e.g., my_bath%Eb(1)%rc%mat) reflect the new values.
    ! This would typically be followed by recalculating the hybridization function:
    ! CALL bath2hybrid(my_bath, FERMIONIC) from bath_class_hybrid module

    CALL delete_bath(my_bath)
  END SUBROUTINE test_bath_vectorization

END MODULE example_bath_vec_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `bath_class`: This module is entirely dependent on the `bath_type` definition from `bath_class.f90`.
    *   `common_def`: Uses `c2s`, `dump_message`, `i2c` for error reporting.
    *   `genvar`: For `DP` (double precision kind).
    *   `globalvar_ed_solver`: Uses the `all_first_call` logical flag to control the `first_call` static variable behavior, ensuring consistency checks run appropriately.
    *   `masked_matrix_class`: Relies heavily on this module, specifically the `masked_matrix_type` and its subroutines `masked_matrix2vec` (to get the 1D vector form of a single matrix's parameters) and `vec2masked_matrix` (to reconstruct a matrix from its 1D vector form).
    *   `masked_matrix_class_mod` (used in `extract_mmvec_from_vec`): For `gather_diag_offdiag_vec`. This suggests a close relationship or potential refactoring opportunity with `masked_matrix_class`. *(Note: The name `masked_matrix_class_mod` seems unusual if `masked_matrix_class` is the primary module; it might indicate a patch or extension module.)*

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit third-party external libraries are directly used.

*   **Interactions with other components:**
    *   This module serves as a critical link between the physically structured representation of bath parameters (`bath_type` in `bath_class`) and numerical optimization routines that require a simple, flat 1D vector of real parameters.
    *   The `bath2vec` routine is called to prepare the parameter vector for an optimizer.
    *   The `vec2bath` routine is called to take an optimized (or trial) parameter vector from an optimizer and update the actual `bath_type` matrices.
    *   It directly calls routines from `masked_matrix_class` to handle the vectorization/reconstruction of individual `masked_matrix_type` objects, correctly managing complex numbers (by splitting into real and imaginary parts) and accounting for hermiticity to ensure only independent parameters are included in the vector.
    *   The module `bath_class_hybrid` would be a primary user of `bath_class_vec` when setting up calls to minimization routines (e.g., `minimize_func_wrapper`).

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
