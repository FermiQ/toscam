# `src/ed_solver/aim2_class.f90`

## Overview

*   **Purpose:** Defines a Fortran derived type `AIM2_type` and associated procedures to represent the quadratic (non-interacting or \(H_0\)) part of an Anderson Impurity Model (AIM). This includes impurity energies, bath energies, and impurity-bath couplings.
*   **Role in Project:** Provides a data structure for the purely quadratic terms of an AIM, effectively separating them from the many-body interaction terms (which would be handled elsewhere). It is designed to handle both normal (spin-polarized) and Nambu (superconducting) representations of the model. This module likely serves as a component within a larger AIM framework or an Exact Diagonalization (ED) solver system.

## Key Components

*   **`TYPE AIM2_type`**:
    A Fortran derived type representing the quadratic part of an Anderson Impurity Model.
    *   `INTEGER :: spin`: Stores the spin sector (e.g., 1 for spin up, 2 for spin down) if the model is not in a superconducting representation. Defaults to 0.
    *   `INTEGER :: Ns`: Total number of sites (impurity + bath).
    *   `INTEGER :: Nc`: Number of impurity sites.
    *   `INTEGER :: Nb`: Number of bath sites.
    *   `LOGICAL :: SUPER`: A boolean flag, `.TRUE.` if the model is in a superconducting (Nambu) representation, `.FALSE.` otherwise.
    *   `TYPE(masked_matrix_type) :: Ec`: A `masked_matrix_type` object (likely defined in `masked_matrix_class`) storing the impurity energy matrix (on-site energies).
    *   `TYPE(masked_matrix_type) :: Eb`: A `masked_matrix_type` object storing the bath energy matrix.
    *   `TYPE(masked_matrix_type) :: Vbc`: A `masked_matrix_type` object storing the coupling (hybridization) matrix between bath and impurity sites.
    *   `INTEGER :: IMPnorbs`: Number of impurity orbitals in the current representation (e.g., `Nc` for normal state, `2*Nc` for Nambu state).
    *   `INTEGER, POINTER :: IMPiorb(:) => NULL()`: A pointer to an integer array that defines the mapping or ordering of the impurity orbitals within the full AIM basis.
    *   `INTEGER :: BATHnorbs`: Number of bath orbitals in the current representation.
    *   `INTEGER, POINTER :: BATHiorb(:) => NULL()`: A pointer to an integer array defining the mapping or ordering of the bath orbitals.

*   **`SUBROUTINE new_AIM2(AIM2, AIM, SPIN)`**:
    The constructor for an `AIM2_type` object. It initializes `AIM2` based on a full `AIM_type` object (`AIM`, presumably from `aim_class.f90`).
    *   **Arguments:**
        *   `AIM2 (TYPE(AIM2_type), INTENT(INOUT))`: The `AIM2_type` object to be initialized.
        *   `AIM (TYPE(AIM_type), INTENT(IN))`: The full Anderson Impurity Model object containing all parameters.
        *   `SPIN (INTEGER, OPTIONAL, INTENT(IN))`: Specifies the spin sector to extract if the model is not superconducting. Required if `AIM%bath%SUPER` is `.FALSE.`.
    *   **Functionality:**
        *   Copies scalar properties (`Nc`, `Nb`, `Ns`, `SUPER`) from `AIM`.
        *   Allocates and defines `IMPiorb` and `BATHiorb` based on the orbital definitions in `AIM`.
        *   If `AIM2%SUPER` is `.TRUE.`, it calls Nambu transformation routines (e.g., `Nambu_Ec`, `Nambu_Eb`, `Nambu_Vbc`, likely imported from `impurity_class` and `bath_class`) to construct the Nambu representation of `Ec`, `Eb`, and `Vbc`.
        *   If `AIM2%SUPER` is `.FALSE.`, it extracts the specified `SPIN` component from the corresponding matrices in `AIM%impurity%Ec(SPIN)`, `AIM%bath%Eb(SPIN)`, and `AIM%bath%Vbc(SPIN)`.
        *   Finally, it updates the masks within the `Ec`, `Eb`, and `Vbc` `masked_matrix_type` objects to reflect only the non-zero matrix elements.

*   **`SUBROUTINE copy_AIM2(AIM2OUT, AIM2IN)`**:
    Performs a deep copy of an `AIM2_type` object from `AIM2IN` to `AIM2OUT`. This includes copying all scalar members and performing deep copies of the `IMPiorb` and `BATHiorb` arrays and the `masked_matrix_type` components (`Ec`, `Eb`, `Vbc`) using their respective copy constructors.

*   **`SUBROUTINE delete_AIM2(AIM2)`**:
    The destructor for an `AIM2_type` object.
    *   Deallocates the `IMPiorb` and `BATHiorb` pointer arrays if they are associated.
    *   Calls the destructor (`delete_masked_matrix`) for each of the `masked_matrix_type` components (`Ec`, `Eb`, `Vbc`) to free their internal memory.

## Important Variables/Constants

Within the `AIM2_type` derived type:

*   **`SUPER :: LOGICAL`**: This is a critical flag that dictates the representation of the model. If `.TRUE.`, matrices are set up in Nambu space for superconducting calculations. If `.FALSE.`, a specific spin sector is typically chosen.
*   **`Ec :: TYPE(masked_matrix_type)`**: Stores the impurity on-site energies. The mask allows for efficient representation of sparse or structured matrices.
*   **`Eb :: TYPE(masked_matrix_type)`**: Stores the bath site energies.
*   **`Vbc :: TYPE(masked_matrix_type)`**: Stores the hybridization (coupling) strengths between impurity and bath sites.
*   **`IMPiorb(:), BATHiorb(:) :: INTEGER, POINTER`**: These arrays define the global indexing of the impurity and bath orbitals, which is crucial for constructing the full Hamiltonian matrix in a consistent manner.

## Usage Examples

```fortran
MODULE example_usage
  USE AIM2_class
  USE AIM_class          ! Assuming AIM_type is defined here
  USE impurity_class     ! For impurity_type
  USE bath_class         ! For bath_type
  USE masked_matrix_class

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: test_aim2

CONTAINS

  SUBROUTINE test_aim2
    TYPE(AIM_type) :: full_aim_model
    TYPE(AIM2_type) :: quadratic_part_spin_up
    TYPE(AIM2_type) :: quadratic_part_nambu

    ! --- Initialize full_aim_model (details omitted) ---
    ! Example: full_aim_model%Nc = 1, full_aim_model%bath%Nb = 2
    ! full_aim_model%Ns = 3
    ! Allocate and set full_aim_model%impurity%Ec, full_aim_model%bath%Eb, etc.
    ! full_aim_model%bath%SUPER = .FALSE. ! For normal state example
    ! CALL new_impurity(...)
    ! CALL new_bath(...)
    ! CALL new_AIM(full_aim_model, ...)
    ! --- End of full_aim_model initialization ---


    ! 1. Create an AIM2 object for the spin-up sector of a normal state AIM
    IF (.NOT. full_aim_model%bath%SUPER) THEN
      CALL new_AIM2(quadratic_part_spin_up, full_aim_model, SPIN=1) ! Assuming spin 1 is up
      PRINT *, "Created AIM2 for normal state, spin up"
      PRINT *, "Impurity orbitals:", quadratic_part_spin_up%IMPnorbs
      ! ... (access quadratic_part_spin_up%Ec%rc%mat, etc.) ...
      CALL delete_AIM2(quadratic_part_spin_up)
    END IF

    ! 2. Create an AIM2 object for a superconducting AIM
    !    (first, conceptually, reconfigure full_aim_model for SUPER = .TRUE.)
    !    full_aim_model%bath%SUPER = .TRUE.
    !    CALL delete_AIM(full_aim_model) ! Clean up previous state
    !    ! ... re-initialize full_aim_model as superconducting ...
    !    CALL new_AIM(full_aim_model, ...)
    !
    !    IF (full_aim_model%bath%SUPER) THEN
    !      CALL new_AIM2(quadratic_part_nambu, full_aim_model) ! SPIN is not given
    !      PRINT *, "Created AIM2 for superconducting state (Nambu)"
    !      PRINT *, "Impurity orbitals (Nambu):", quadratic_part_nambu%IMPnorbs
    !      ! ... (access quadratic_part_nambu%Ec%rc%mat - now a Nambu matrix) ...
    !      CALL delete_AIM2(quadratic_part_nambu)
    !    END IF

    CALL delete_AIM(full_aim_model) ! Clean up full model

  END SUBROUTINE test_aim2

END MODULE example_usage
```
This module is primarily intended for internal use by other modules within an ED or DMFT solver that need to work with the non-interacting part of an Anderson Impurity Model.

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `masked_matrix_class`: Crucially depends on this module for the `masked_matrix_type` definition and its associated procedures (`new_masked_matrix`, `copy_masked_matrix`, `delete_masked_matrix`). `AIM2_type` uses `masked_matrix_type` for `Ec`, `Eb`, and `Vbc`.
    *   `aim_class` (specifically in `new_AIM2`): The constructor `new_AIM2` takes an object of `AIM_type` (defined in `aim_class.f90`) as input to extract the quadratic Hamiltonian parameters.
    *   `bath_class` (specifically in `new_AIM2`): Uses helper functions `nambu_eb` and `nambu_vbc` from this module when constructing the Nambu representation for a superconducting AIM.
    *   `impurity_class` (specifically in `new_AIM2`): Uses the helper function `nambu_ec` from this module for the Nambu representation of impurity energies.
    *   `genvar` (specifically in `new_AIM2`): Uses the `DP` constant for defining double precision real kinds.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit external (third-party) libraries are directly linked or used by this module.

*   **Interactions with other components:**
    *   This module acts as a data structure and utility provider. It is expected to be used by higher-level modules within an ED or DMFT solver framework.
    *   The `new_AIM2` subroutine is the primary point of interaction with other data structures like `AIM_type`. It effectively "distills" the quadratic, non-interacting part from a more complete impurity model description.
    *   The resulting `AIM2_type` object can then be used by other routines to, for example:
        *   Construct the matrix representation of the non-interacting part of the Hamiltonian.
        *   Calculate non-interacting Green's functions.
        *   Serve as input for routines that add interaction terms to build the full Hamiltonian.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
