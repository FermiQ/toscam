# `src/ed_solver/aim_class.f90`

## Overview

*   **Purpose:** Defines a Fortran derived type `AIM_type` and associated procedures to represent a complete Anderson Impurity Model (AIM). It encapsulates both impurity and bath characteristics, defines their combined structure, and provides utilities for managing state representations within this combined system.
*   **Role in Project:** This module serves as a central data structure for defining an Anderson Impurity Model. It acts as a high-level container that links an `impurity_type` object and a `bath_type` object. It establishes a consistent framework for the total number of sites, orbitals, and the overall Hilbert space size. Furthermore, it includes crucial helper routines to map integer indices between the combined AIM many-body states and the individual impurity and bath states, which is essential for constructing and interpreting these states in Exact Diagonalization (ED) calculations.

## Key Components

*   **`TYPE AIM_type`**:
    A Fortran derived type representing the full Anderson Impurity Model.
    *   `TYPE(impurity_type), POINTER :: impurity => NULL()`: A pointer to an `impurity_type` object, which contains detailed information about the impurity sites (e.g., energies, interactions).
    *   `TYPE(bath_type), POINTER :: bath => NULL()`: A pointer to a `bath_type` object, which describes the non-interacting bath coupled to the impurity (e.g., bath site energies, hybridization values).
    *   `INTEGER :: Ns`: Total number of sites in the AIM (calculated as `Nc + Nb`).
    *   `INTEGER :: norbs`: Total number of single-particle orbitals in the AIM representation (e.g., for a spinful model, this would typically be `(Nc + Nb) * 2`).
    *   `INTEGER :: nstates`: Total number of many-body states in the combined Hilbert space. Initialized as `impurity%nstates * bath%nstates`, suggesting it represents the dimension of the product space of impurity and bath states.
    *   `INTEGER :: Nc`: Number of impurity sites, derived from `impurity%Nc`.
    *   `INTEGER :: Nb`: Number of bath sites, derived from `bath%Nb`.
    *   `INTEGER, POINTER :: IMPiorb(:, :) => NULL()`: A 2D allocatable array (`site_index, spin_index`) that stores the global index of each impurity orbital within the combined AIM single-particle basis.
    *   `INTEGER, POINTER :: BATHiorb(:, :) => NULL()`: A 2D allocatable array (`site_index, spin_index`) that stores the global index of each bath orbital within the combined AIM single-particle basis.

*   **`SUBROUTINE new_AIM(AIM, impurity, bath)`**:
    The constructor for an `AIM_type` object.
    *   **Arguments:**
        *   `AIM (TYPE(AIM_type), INTENT(INOUT))`: The `AIM_type` object to be initialized.
        *   `impurity (TYPE(impurity_type), INTENT(IN), TARGET)`: The target `impurity_type` object.
        *   `bath (TYPE(bath_type), INTENT(IN), TARGET)`: The target `bath_type` object.
    *   **Functionality:**
        *   Calls `delete_AIM` to ensure a clean state.
        *   Calls `update_AIM_pointer` to associate `AIM%impurity` and `AIM%bath` with the provided `impurity` and `bath` objects.
        *   Sets scalar properties: `Nc`, `Nb`, `Ns`, `norbs`, `nstates`.
        *   Allocates `IMPiorb` and `BATHiorb` (both as `(max_sites_in_type, 2)` arrays, e.g., `(AIM%Nb, 2)` for bath).
        *   Initializes these orbital mapping arrays. The documented ordering convention is:
            1.  Bath orbitals, spin up
            2.  Impurity orbitals, spin up
            3.  Bath orbitals, spin down
            4.  Impurity orbitals, spin down
            The `ramp` function (from `linalg`) is used to assign sequential indices based on this order.

*   **`SUBROUTINE update_AIM_pointer(AIM, impurity, bath)`**:
    A utility subroutine that associates the `AIM%impurity` and `AIM%bath` pointers with the respective target `impurity_type` and `bath_type` objects.

*   **`SUBROUTINE delete_AIM(AIM)`**:
    The destructor for an `AIM_type` object.
    *   Deallocates the `BATHiorb` and `IMPiorb` arrays if they are associated (using `istati` from `globalvar_ed_solver` for status).
    *   Nullifies the `AIM%impurity` and `AIM%bath` pointers.

*   **`SUBROUTINE IMPBATH2AIMstate(AIMstate, IMPstate, BATHstate, Nc, Nb)`**:
    Converts separate integer indices representing an impurity state (`IMPstate`) and a bath state (`BATHstate`) into a single integer index for the corresponding combined AIM state (`AIMstate`).
    *   **Logic:** It first resolves the up and down spin components of the `BATHstate` and `IMPstate` based on the number of states in their respective single-spin subspaces (`nBATHstates_sz = 2**Nb`, `nIMPstates_sz = 2**Nc`). Then, it combines these components according to the defined orbital ordering (bath-up, imp-up, then bath-down, imp-down) to form the `AIMstate` index.

*   **`SUBROUTINE AIM2IMPBATHstate(IMPstate, BATHstate, AIMstate, Nc, Nb)`**:
    The inverse operation of `IMPBATH2AIMstate`. It takes a combined `AIMstate` index and decomposes it into the corresponding `IMPstate` and `BATHstate` indices.
    *   **Logic:** It reverses the steps of `IMPBATH2AIMstate`, first separating the global `AIMstate` into its up and down spin components, then further decomposing each of these into their bath and impurity parts.

## Important Variables/Constants

Within the `AIM_type` derived type:

*   **`impurity :: TYPE(impurity_type), POINTER`**: A pointer to the `impurity_type` object, which contains all details of the impurity sites, local interactions, etc.
*   **`bath :: TYPE(bath_type), POINTER`**: A pointer to the `bath_type` object, which describes the non-interacting bath sites and their coupling to the impurity.
*   **`Ns, Nc, Nb, norbs, nstates :: INTEGER`**: These scalar integers define the overall size and complexity of the Anderson Impurity Model's Hilbert space.
*   **`IMPiorb(:,:), BATHiorb(:,:) :: INTEGER, POINTER`**: These arrays are critical for defining the specific mapping of local impurity and bath orbitals to a global indexing scheme for the combined AIM. This ordering is essential for correctly constructing Hamiltonian matrices and applying operators in the combined basis.

## Usage Examples

```fortran
MODULE example_aim_usage
  USE AIM_class
  USE impurity_class, ONLY: impurity_type, new_impurity, delete_impurity ! Assuming these exist
  USE bath_class, ONLY: bath_type, new_bath, delete_bath       ! Assuming these exist
  USE linalg, ONLY: ramp ! For initializing some parts of impurity/bath if needed

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_aim_module

CONTAINS

  SUBROUTINE test_aim_module
    TYPE(impurity_type), TARGET :: test_impurity
    TYPE(bath_type), TARGET :: test_bath
    TYPE(AIM_type) :: my_full_aim

    INTEGER :: imp_s, bath_s, aim_s_combined
    INTEGER :: imp_s_test, bath_s_test

    ! Initialize test_impurity and test_bath (simplified example)
    test_impurity%Nc = 1
    test_impurity%norbs = 2 * test_impurity%Nc ! Spinful
    test_impurity%nstates = 2**test_impurity%norbs ! Incorrect simplified nstates, just for demo
                                                 ! Actual nstates is usually 2**Nc for impurity part in product basis
    test_bath%Nb = 2
    test_bath%norbs = 2 * test_bath%Nb ! Spinful
    test_bath%nstates = 2**test_bath%norbs ! Incorrect simplified nstates

    ! In a real scenario, new_impurity and new_bath would populate these properly
    ! CALL new_impurity(test_impurity, ...)
    ! CALL new_bath(test_bath, ...)

    ! Create the AIM object
    CALL new_AIM(my_full_aim, test_impurity, test_bath)

    PRINT *, "AIM Model Created:"
    PRINT *, "  Total sites (Ns): ", my_full_aim%Ns
    PRINT *, "  Impurity sites (Nc): ", my_full_aim%Nc
    PRINT *, "  Bath sites (Nb): ", my_full_aim%Nb
    PRINT *, "  Total orbitals (norbs): ", my_full_aim%norbs
    PRINT *, "  Total states (nstates from product): ", my_full_aim%nstates
    PRINT *, "  Impurity orbital mapping (1st site, up/down): ", my_full_aim%IMPiorb(1,:)
    PRINT *, "  Bath orbital mapping (1st site, up/down): ", my_full_aim%BATHiorb(1,:)

    ! Example state conversion (indices are 0-based in example, adjust if 1-based)
    imp_s = 1  ! Example impurity state index (0 to 2**Nc - 1)
    bath_s = 2 ! Example bath state index (0 to 2**Nb - 1)

    CALL IMPBATH2AIMstate(aim_s_combined, imp_s, bath_s, my_full_aim%Nc, my_full_aim%Nb)
    PRINT *, "For IMPstate=", imp_s, " and BATHstate=", bath_s, ", combined AIMstate=", aim_s_combined

    CALL AIM2IMPBATHstate(imp_s_test, bath_s_test, aim_s_combined, my_full_aim%Nc, my_full_aim%Nb)
    PRINT *, "Converted back: IMPstate=", imp_s_test, ", BATHstate=", bath_s_test

    ! Clean up
    CALL delete_AIM(my_full_aim)
    ! CALL delete_impurity(test_impurity) ! If dynamically allocated
    ! CALL delete_bath(test_bath)       ! If dynamically allocated

  END SUBROUTINE test_aim_module

END MODULE example_aim_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `bath_class`: Essential for the `bath_type` definition, which is a component of `AIM_type`.
    *   `impurity_class`: Essential for the `impurity_type` definition, another component of `AIM_type`.
    *   `linalg` (used in `new_AIM`): The `ramp` subroutine is used for initializing orbital mapping arrays (`IMPiorb`, `BATHiorb`).
    *   `globalvar_ed_solver` (used in `delete_AIM`): The `istati` variable is likely an integer status variable used during deallocation to check for errors.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit third-party external libraries are directly used by this module.

*   **Interactions with other components:**
    *   This module serves as a high-level container that aggregates `impurity_type` and `bath_type` objects. It is fundamental for ED solvers and other routines that need a complete description of the AIM.
    *   The orbital mapping arrays `IMPiorb` and `BATHiorb` are critical for any routine that needs to construct Hamiltonian matrices or apply operators in the combined impurity+bath basis, ensuring consistent indexing.
    *   The state conversion routines (`IMPBATH2AIMstate`, `AIM2IMPBATHstate`) are particularly important for ED algorithms that work by constructing many-body basis states as products of impurity states and bath states. These routines allow for seamless translation between the "local" state representations and the "global" AIM state representation.
    *   Other modules, such as `AIM2_class` (which defines the quadratic part of an AIM), would typically take an `AIM_type` object as input to extract the information they need.
    *   ED solver routines would use the `AIM_type` to understand the structure of the problem, build Hamiltonian matrices, and interpret the resulting eigenstates.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
