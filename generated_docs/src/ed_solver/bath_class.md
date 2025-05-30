# `src/ed_solver/bath_class.f90`

## Overview

*   **Purpose:** Defines a Fortran derived type `bath_type` and associated procedures to represent the non-interacting electronic bath that is coupled to a quantum impurity within the framework of an Anderson Impurity Model (AIM).
*   **Role in Project:** This module is responsible for managing all aspects of the bath component of an AIM. This includes storing bath site energies (`Eb`), the hybridization strengths coupling bath sites to impurity sites (`Vbc`), and, for superconducting scenarios, pairing terms within the bath (`Pb`) and pairing hybridization (`PVbc`). It provides a comprehensive suite of functionalities for creating, initializing (from scratch or files), copying, deleting, reading, and writing these bath definitions. Furthermore, it includes routines to transform the bath representation into Nambu space for superconducting calculations and stores parameters related to the process of fitting the bath model to a target hybridization function.

## Key Components

*   **`TYPE bath_type`**:
    A Fortran derived type encapsulating the properties of the electronic bath.
    *   `CHARACTER(LEN=100) :: fileout`: Default filename for writing bath parameters.
    *   `INTEGER :: Nb`: Number of bath sites.
    *   `INTEGER :: Nc`: Number of impurity sites (relevant for defining hybridization `Vbc`).
    *   `INTEGER :: norbs`: Total number of single-particle orbitals in the bath (typically `Nb*2` for a spinful bath).
    *   `INTEGER :: nstates`: Dimension of the bath's Hilbert space (typically `2**norbs`).
    *   `INTEGER, POINTER :: iorb(:, :) => NULL()`: A 2D array mapping `(site_index, spin_index)` for bath orbitals to a global index.
    *   `LOGICAL :: SUPER = .FALSE.`: A flag indicating if the bath is in a superconducting (Nambu) representation. Determined by `which_basis_to_use`.
    *   `TYPE(masked_matrix_type), POINTER :: Eb(:) => NULL()`: An array (size 2, for spin-up and spin-down) of `masked_matrix_type` objects storing bath site energies.
    *   `TYPE(masked_matrix_type), POINTER :: Vbc(:) => NULL()`: An array (size 2) for the hybridization matrices \(V_{bath,impurity}\) for each spin.
    *   `TYPE(masked_matrix_type), POINTER :: PVbc(:) => NULL()`: An array (size 2) for pairing hybridization matrices (e.g., \(\Delta_{bath,impurity}\)) coupling bath and impurity, relevant for superconducting states.
    *   `TYPE(masked_matrix_type) :: Pb`: A `masked_matrix_type` object for pairing terms *within* the bath itself (e.g., \(\Delta_{bath,bath}\)).
    *   `INTEGER :: nparam`: Total number of independent real parameters defining the bath Hamiltonian terms.
    *   `REAL(DP), POINTER :: vec(:) => NULL()`: A 1D array that can store all bath parameters in a flat vector, potentially used for optimization/fitting.
    *   `INTEGER :: Niter_search_max`: Maximum number of iterations for fitting the bath parameters to a target hybridization function.
    *   `REAL(DP) :: search_step`: Initial step size used in the hybridization fitting search algorithm.
    *   `REAL(DP) :: dist_max`: Maximum allowed error or distance when fitting the hybridization function.
    *   `TYPE(correl_type) :: hybrid, hybridret`: `correl_type` objects (from `correl_class`) storing the Matsubara and retarded hybridization functions \(\Delta(\omega)\) of the bath.

*   **`INTERFACE new_bath`**: A generic interface for creating `bath_type` objects.
    *   `MODULE PROCEDURE new_bath_from_scratch(BATH, Nb, Nc, IMASKE, IMASKP, IMASKV, IMASKPV)`: Creates a new bath object from scratch, defining its structure (number of sites, masks for non-zero matrix elements). `IMASKE`, `IMASKP`, `IMASKV`, `IMASKPV` are optional integer arrays specifying these masks.
    *   `MODULE PROCEDURE new_bath_from_old(BATHOUT, BATHIN)`: Creates `BATHOUT` by first setting up its structure based on `BATHIN` (using `new_bath_from_scratch` with masks from `BATHIN`) and then copying the content from `BATHIN` using `copy_bath`.

*   **`SUBROUTINE copy_bath(BATHOUT, BATHIN)`**: Performs a deep copy of all components from `BATHIN` to `BATHOUT`.
*   **`SUBROUTINE delete_bath(BATH)`**: Destructor for `bath_type` objects, deallocating all dynamically allocated arrays and `masked_matrix_type` components.
*   **`SUBROUTINE read_bath(BATH, FILEIN, Nc)`**: Reads bath parameters from a formatted input file `FILEIN`. This is a complex routine that initializes masks and parameter values based on the input file content and various global flags (e.g., `min_all_bath_param` for random initialization, `diag_bath`, `force_pairing_from_mask`, `para_state`).
*   **`SUBROUTINE read_raw_bath(BATH, FILEIN)`**: Reads bath parameters from a file `FILEIN` assumed to be in a simpler "raw" format.
*   **`SUBROUTINE write_raw_bath(BATH)`**: Writes the current bath parameters to a "raw" format file specified by `BATH%fileout`.
*   **`SUBROUTINE write_bath(BATH, UNIT, FILEIN)`**: Writes a formatted description of the bath parameters to a specified Fortran unit `UNIT` or to the `log_unit`. Can also note the input file `FILEIN`.
*   **`FUNCTION nbathparam(BATH) AS INTEGER`**: Calculates the number of independent real parameters in the bath definition. It accounts for hermiticity and whether matrices are complex or real (via the internal `offset__` function).
*   **`FUNCTION nbathparam_(BATH) AS INTEGER`**: A simpler version to calculate the number of bath parameters by summing the number of non-zero elements specified by the masks (`MASK%nind`).
*   **`SUBROUTINE Nambu_Eb(EbNambu, Eb, Pb)`**: Constructs the Nambu representation of bath energies (`EbNambu`). The diagonal blocks are \(E_{b,\uparrow}\) and \(-E_{b,\downarrow}^T\), and off-diagonal blocks are \(P_b\) and \(P_b^\dagger\).
*   **`SUBROUTINE Nambu_Vbc(VbcNambu, Vbc, PVbc)`**: Constructs the Nambu representation of the hybridization matrix (`VbcNambu`). Diagonal blocks are \(V_{bc,\uparrow}\) and \(-V_{bc,\downarrow}^*\) (or \(-V_{bc,\downarrow}\) if real). Off-diagonal blocks are formed from `PVbc` if `PAIRING_IMP_TO_BATH` is true.
*   **`SUBROUTINE which_basis_to_use(BATH)` / `LOGICAL FUNCTION which_basis_to_use__(mat)`**: Determines if a superconducting basis (`BATH%SUPER = .TRUE.`) should be used. This decision is based on whether pairing terms `Pb` are present in the mask, or if global flags like `supersc_state` or `force_sz_basis` are set true (unless overridden by `force_nupdn_basis`).

## Important Variables/Constants

Within the `bath_type` derived type:

*   **`SUPER :: LOGICAL`**: A critical flag that determines if the bath is treated in a normal state representation or a superconducting (Nambu) representation. This affects how matrices are constructed and interpreted.
*   **`Eb(:), Vbc(:), Pb, PVbc(:) :: TYPE(masked_matrix_type)`**: These `masked_matrix_type` objects store the core Hamiltonian parameters of the bath: `Eb` for on-site bath energies (per spin), `Vbc` for impurity-bath hybridization (per spin), `Pb` for pairing terms within the bath, and `PVbc` for pairing terms between impurity and bath.
*   **`hybrid, hybridret :: TYPE(correl_type)`**: These store the Matsubara and retarded hybridization functions, respectively. These functions are central to DMFT and other embedding methods as they quantify the influence of the bath on the impurity.
*   **`nparam :: INTEGER`**: Stores the total number of adjustable parameters that define the bath. This is used, for example, in fitting procedures where these parameters are optimized.
*   **`Niter_search_max, search_step, dist_max`**: Parameters controlling the numerical procedure for fitting the bath model to a target hybridization function.

## Usage Examples

```fortran
MODULE example_bath_usage
  USE bath_class
  USE correl_class, ONLY: correl_type ! For hybrid definition
  USE masked_matrix_class, ONLY: masked_matrix_type ! For matrix definitions

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_bath_module

CONTAINS
  SUBROUTINE test_bath_module
    TYPE(bath_type) :: my_bath
    INTEGER :: num_imp_sites
    CHARACTER(LEN=100) :: bath_input_file

    num_imp_sites = 1
    bath_input_file = "bath_config.dat" ! This file must exist and be formatted for read_bath

    ! Initialize the bath by reading from a configuration file
    CALL read_bath(my_bath, bath_input_file, num_imp_sites)

    PRINT *, "Bath Initialized. Number of Bath Sites (Nb): ", my_bath%Nb
    PRINT *, "Is superconducting (SUPER): ", my_bath%SUPER
    PRINT *, "Number of parameters (nparam): ", my_bath%nparam

    ! Print bath details to standard output
    CALL write_bath(my_bath)

    ! If the bath is superconducting, one might want to get Nambu matrices
    IF (my_bath%SUPER) THEN
      TYPE(masked_matrix_type) :: nambu_bath_energies
      TYPE(masked_matrix_type) :: nambu_hybridization_coupling

      CALL Nambu_Eb(nambu_bath_energies, my_bath%Eb, my_bath%Pb)
      CALL Nambu_Vbc(nambu_hybridization_coupling, my_bath%Vbc, my_bath%PVbc)

      PRINT *, "Nambu matrices constructed."
      ! (Code to use or inspect nambu_bath_energies and nambu_hybridization_coupling)

      CALL delete_masked_matrix(nambu_bath_energies)
      CALL delete_masked_matrix(nambu_hybridization_coupling)
    END IF

    ! Clean up
    CALL delete_bath(my_bath)

  END SUBROUTINE test_bath_module
END MODULE example_bath_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `correl_class`: For `correl_type`, used to store hybridization functions (`hybrid`, `hybridret`).
    *   `genvar`: Provides global constants like `DP` (double precision kind), `cspin` (spin character), `log_unit`, `rank`.
    *   `masked_matrix_class`: Essential for `masked_matrix_type`, which is used to store all Hamiltonian terms (`Eb`, `Vbc`, `Pb`, `PVbc`).
    *   `globalvar_ed_solver`: Supplies various global flags and default parameters that influence bath initialization and behavior (e.g., `force_nupdn_basis`, `energy_global_shift2`, `pairing_imp_to_bath`, `min_all_bath_param`, `Niter_search_max`, `search_step`, `dist_max`).
    *   `common_def`: For file I/O utilities (`open_safe`, `close_safe`, `skip_line`) and messaging (`dump_message`).
    *   `linalg`: The `ramp` function is used in `new_bath_from_scratch` for initializing orbital indices.
    *   `matrix`: The `write_array` and `diag` functions are used, for example, in `Nambu_Eb`.
    *   `random`: The `dran_tab` function is used in `read_bath` for random initialization of parameters if `min_all_bath_param` is set appropriately.
    *   `string5`: For `get_unit` utility in file operations.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit third-party external libraries are directly linked by this module.

*   **Interactions with other components:**
    *   The `bath_type` object is a fundamental component of the `AIM_type` (defined in `aim_class.f90`), representing the non-interacting bath part of the Anderson Impurity Model.
    *   It relies heavily on `masked_matrix_class` for storing and managing the matrices that define the bath Hamiltonian and its coupling to the impurity.
    *   The `read_bath` subroutine demonstrates complex interaction with global settings (from `globalvar_ed_solver`) to allow flexible and context-dependent initialization of bath parameters.
    *   The Nambu transformation routines (`Nambu_Eb`, `Nambu_Vbc`) are critical for ED solvers that need to work in a superconducting basis.
    *   The hybridization functions (`hybrid`, `hybridret`) stored within `bath_type` are key outputs of bath fitting procedures (likely handled by or in conjunction with `bath_class_hybrid`) and serve as input to the impurity solver.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
