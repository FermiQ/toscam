# `src/ed_solver/impurity_class.f90`

## Overview

*   **Purpose:** This module defines the `impurity_type` Fortran derived type, along with several auxiliary (though not publicly exported from this module) types: `symmetry`, `unitcell`, `web`, and `Hamiltonian`. Together, these describe the quantum impurity within an Anderson Impurity Model (AIM), including its site structure, on-site energies, local Coulomb interactions, and potentially its embedding within a larger lattice context.
*   **Role in Project:** It serves as a central data structure for defining the local interacting region (the "impurity") in an Exact Diagonalization (ED) or Dynamical Mean-Field Theory (DMFT) calculation. The module provides procedures to initialize, define, manage, and write out the properties of this impurity, which are then used by other parts of the solver to construct and solve the AIM.

## Key Components

*   **Auxiliary Derived Types (defined in the module, but not public):**
    *   **`TYPE symmetry`**: Intended to store information about point group symmetries of a crystal or molecule. Includes arrays for rotation matrices (`rotmat`), translation vectors (`trans`, `ttrans`), character tables (`rotchar`), etc.
    *   **`TYPE unitcell`**: Designed to describe the geometry of a unit cell, including lattice vectors (`a`, `b`, `c`), reciprocal lattice vectors, site positions within the cell (`insidepos`), neighbor connectivity (`nneigh`), k-space information (`matkxky`, `dispersion`), etc.
    *   **`TYPE web`**: Combines `unitcell` and `symmetry` types. Appears to represent a larger, potentially periodic lattice structure from which the impurity sites can be selected. Contains information about links between sites (`struct_link`), site positions (`struct_pos`), and various geometric and topological properties.
    *   **`TYPE Hamiltonian`**: A general structure for holding various terms of a Hamiltonian, such as hopping parameters (`teta`), pairing terms (`delta`), on-site energies (`eps`), Coulomb repulsion terms (`Vrep`), Hund's coupling (`hund`), and external fields (`field`).

*   **`TYPE impurity_type` (Public):**
    The primary data structure representing the quantum impurity.
    *   `INTEGER :: Nc = 0`: The number of sites that constitute the impurity.
    *   `INTEGER :: norbs = 0`: The total number of single-particle orbitals on the impurity (e.g., if each site is spinful, `norbs = Nc * 2`).
    *   `INTEGER :: nstates = 0`: The dimension of the local Hilbert space for the impurity (typically \(2^{\text{norbs}}\)).
    *   `INTEGER, POINTER :: iorb(:, :) => NULL()`: A 2D array, typically `(Nc, 2)`, that maps a local impurity site index and a spin index to a global or internal orbital index scheme.
    *   `TYPE(masked_matrix_type), POINTER :: Ec(:) => NULL()`: An allocatable array (usually size 2, for spin-up and spin-down components) of `masked_matrix_type` objects. These store the quadratic part of the impurity Hamiltonian, including on-site energies and local magnetic fields.
    *   `TYPE(masked_real_matrix_type) :: U`: A `masked_real_matrix_type` object storing the local Coulomb interaction terms on the impurity (e.g., Hubbard U values on each site, inter-site V values between impurity sites). This matrix is real and symmetric.

*   **Public Module Variable:**
    *   `TYPE(masked_matrix_type), PUBLIC, SAVE :: Eccc`: A saved `masked_matrix_type` variable, which is used to store the Nambu representation of the impurity's quadratic Hamiltonian terms after a call to `update_impurity` (which internally calls `Nambu_Ec`).

*   **Public Subroutines:**
    *   **`SUBROUTINE new_impurity(impurity, Nc, IMASKE, IMASKU)`**: The constructor for `impurity_type`.
        *   Initializes `impurity` for `Nc` sites. Sets `norbs = Nc*2` and `nstates = 2**norbs`.
        *   Initializes `iorb` to map site `i` spin up to `i` and spin down to `i+Nc`.
        *   Optionally, if `IMASKE` (for `Ec`) and `IMASKU` (for `U`) integer masks are provided, it initializes the `Ec` and `U` `masked_matrix_type` structures with these masks.
    *   **`SUBROUTINE delete_impurity(IMP)`**: Destructor for `impurity_type`. Deallocates `iorb` and the `Ec` array (calling `delete_masked_matrix` for each element). The `U` matrix is not explicitly deallocated here, assuming its `masked_real_matrix_type` handles its own memory or is deallocated elsewhere if it's a pointer.
    *   **`SUBROUTINE copy_impurity(IMPOUT, IMPIN)`**: Performs a deep copy of an `impurity_type` from `IMPIN` to `IMPOUT`.
    *   **`SUBROUTINE define_impurity(impurity, mmu, impurity_, Himp, Eimp)`**: Populates the `impurity_type` object with more detailed physical parameters.
        *   It first calls `delete_impurity` and then `new_impurity`.
        *   The Coulomb interaction matrix `impurity%U` is set up using parameters from `Himp%dU` (on-site part) and `Himp%Vrep` (inter-site part, if `impurity_%cadran` indicates such links), or from the global `UUmatrix` if available.
        *   It then calls `update_impurity` to set up the quadratic terms `impurity%Ec`.
    *   **`SUBROUTINE update_impurity(Nc, impurity, mmu, impurity_, Himp, Eimp)`**: Updates the quadratic part (`Ec`) of the `impurity_type`.
        *   If `Eimp` (an explicit impurity energy matrix) is provided, `Ec` is derived from `Eimp` (with chemical potential `mmu` subtracted from diagonal elements). The Nambu form for spin-down is \( -(E_{imp, \downarrow\downarrow})^T \).
        *   If `Eimp` is not provided, `Ec` is derived from `Himp%eps` (on-site energies minus `mmu`) and `Himp%teta` (local hoppings/fields within the impurity sites defined by `impurity_`).
        *   After `Ec` is set up, it calls `Nambu_Ec` to compute and store the Nambu version in the public module variable `Eccc`.
    *   **`SUBROUTINE write_impurity(impurity, UNIT)`**: Prints a formatted summary of the impurity parameters (`Nc`, `Ec` matrices, `U` matrix) to a specified Fortran `UNIT` or the default log unit.
    *   **`SUBROUTINE Nambu_Ec(EcNambu, Ec)`**: Constructs the Nambu representation `EcNambu` (a single `masked_matrix_type`) from the spin-dependent `Ec(:)` array. The top-left block is `Ec(1)%rc%mat` (spin up), and the bottom-right block is `-TRANSPOSE(Ec(2)%rc%mat)` (related to spin down). It also calculates and stores a global energy shift (`energy_global_shift` from `globalvar_ed_solver`).

*   **Other Utility Functions (not all public from the module scope, but defined within):**
    *   `FUNCTION average_chem_pot(impurity) AS REAL(DP)`: Calculates the average chemical potential from the diagonal elements of `impurity%Ec`.
    *   `SUBROUTINE shift_average_chem_pot(new_chem_pot, impurity)`: Shifts the diagonal elements of `impurity%Ec` such that their average becomes `new_chem_pot`.
    *   `SUBROUTINE T1_T2_connect_unitcells(...)`: A helper for the `web` type, determining unit cell connectivity based on lattice vectors.

## Important Variables/Constants

*   **Within `impurity_type`:**
    *   `Nc :: INTEGER`: The number of physical sites comprising the impurity.
    *   `Ec(:) :: TYPE(masked_matrix_type), POINTER`: Array (typically size 2 for spin-up, spin-down) storing the quadratic Hamiltonian terms for the impurity (on-site energies, local magnetic fields, intra-impurity hoppings).
    *   `U :: TYPE(masked_real_matrix_type)`: Stores the local Coulomb interaction parameters (e.g., Hubbard U on each site, inter-site density-density interactions V within the impurity cluster).
*   **Module Variable `Eccc :: TYPE(masked_matrix_type), PUBLIC, SAVE`**: Stores the Nambu representation of the quadratic impurity Hamiltonian terms, updated whenever `update_impurity` is called.
*   **Auxiliary types `web` and `Hamiltonian`**: While not public, these are used by `define_impurity` to transfer information from a more general lattice and Hamiltonian description into the specific `impurity_type` structure.

## Usage Examples

```fortran
MODULE example_impurity_usage
  USE impurity_class
  USE masked_matrix_class, ONLY: delete_masked_matrix ! For Eccc if used
  USE genvar, ONLY: DP

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: setup_and_inspect_impurity

CONTAINS
  SUBROUTINE setup_and_inspect_impurity
    TYPE(impurity_type) :: my_imp_problem
    INTEGER, PARAMETER :: num_impurity_sites = 1
    REAL(DP) :: chemical_potential

    ! --- Option 1: Basic initialization ---
    CALL new_impurity(my_imp_problem, Nc=num_impurity_sites)
    ! At this point, my_imp_problem%Ec and my_imp_problem%U are allocated
    ! (if IMASKE/IMASKU were passed or if new_masked_matrix defaults to allocating).
    ! Their values would need to be filled, e.g.:
    ! IF (ASSOCIATED(my_imp_problem%Ec(1)%rc%mat)) my_imp_problem%Ec(1)%rc%mat(1,1) = -0.5_DP ! Example on-site energy
    ! IF (ASSOCIATED(my_imp_problem%U%mat)) my_imp_problem%U%mat(1,1) = 4.0_DP           ! Example Hubbard U

    ! --- Option 2: Define from more general 'web' and 'Hamiltonian' types (conceptual) ---
    ! TYPE(web) :: lattice_subset_for_impurity
    ! TYPE(Hamiltonian) :: system_hamiltonian
    ! chemical_potential = 0.0_DP
    ! ! (Code to initialize lattice_subset_for_impurity and system_hamiltonian)
    ! CALL define_impurity(my_imp_problem, chemical_potential, lattice_subset_for_impurity, system_hamiltonian)
    ! This internally calls update_impurity which also sets the module variable Eccc.

    ! --- Or directly update Ec if U is already set (e.g. after new_impurity) ---
    ! chemical_potential = 0.0_DP
    ! TYPE(web) :: dummy_web; dummy_web%N = num_impurity_sites ! Simplified for example
    ! TYPE(Hamiltonian) :: dummy_H; ! Allocate H%eps if needed
    ! CALL update_impurity(num_impurity_sites, my_imp_problem, chemical_potential, dummy_web, dummy_H)


    PRINT *, "Impurity Details:"
    PRINT *, "  Number of sites (Nc): ", my_imp_problem%Nc
    CALL write_impurity(my_imp_problem)

    ! Accessing the Nambu form of Ec (Eccc is a public module variable)
    IF (ASSOCIATED(Eccc%rc%mat)) THEN
      PRINT *, "Nambu form of Ec (Eccc) is available."
      ! CALL write_masked_matrix(Eccc)
    END IF

    ! Clean up
    CALL delete_impurity(my_imp_problem)
    IF (ASSOCIATED(Eccc%rc%mat)) THEN
        CALL delete_masked_matrix(Eccc) ! Important to clean up saved module variable if appropriate
    END IF

  END SUBROUTINE setup_and_inspect_impurity
END MODULE example_impurity_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP` (double precision kind), `cspin` (character spin labels), `log_unit`, `rank`.
    *   `globalvar_ed_solver`: For `istati` (status variable for allocations) and global parameters like `UUmatrix` (default U values) and `energy_global_shift` (modified by `Nambu_Ec`).
    *   `linalg`: For `ramp` (used in `new_impurity` to initialize `iorb`) and `arrayi` (type used in `symmetry`). The `decomposevec` and `cells` routines are used by the private `T1_T2_connect_unitcells`.
    *   `masked_matrix_class`: This is a crucial dependency for `masked_matrix_type` (used for `Ec`) and `masked_real_matrix_type` (used for `U`), along with their associated procedures for creation, deletion, copying, filling, and testing.
    *   `matrix` (used in `define_impurity`, `Nambu_Ec`, `average_chem_pot`): For `write_array` (printing matrices) and `diag` (extracting diagonal elements).

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.

*   **Interactions with other components:**
    *   **`AIM_class.f90`**: The `impurity_type` object is a fundamental component of the higher-level `AIM_type` (Anderson Impurity Model type). `AIM_type` typically holds a pointer to an `impurity_type` object.
    *   **Hamiltonian Construction (e.g., `h_class.f90`):** Modules responsible for constructing the full Hamiltonian matrix for the AIM will consume `impurity_type` objects to get the on-site energies (`Ec`) and local interaction terms (`U`).
    *   **ED Solvers (e.g., `dmft_solver_ed.f90`):** The main solver routines will directly or indirectly use the `impurity_type` to define the problem being solved. `update_impurity` is often called within solver loops if parameters like chemical potential or explicit impurity energies `Eimp` change.
    *   **Nambu Formalism:** The `Nambu_Ec` routine and the public module variable `Eccc` are key for calculations involving superconductivity, where the Hamiltonian needs to be expressed in particle-hole (Nambu) space.
    *   The auxiliary types (`web`, `Hamiltonian`) suggest that this module can interface with codes that define impurities based on more complex lattice geometries or general Hamiltonian parameterizations, providing a pathway to extract the specific impurity information needed by the ED solver.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
