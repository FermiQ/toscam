# `src/ed_solver/apply_c.f90`

## Overview

*   **Purpose:** This module provides a suite of Fortran subroutines designed to apply fermionic creation (\(c^\dagger\)) or destruction (\(c\)) operators for a specified spin (up/down) and impurity site to a given many-body eigenstate. Such operations are fundamental for calculating Green's functions, susceptibilities, and other multi-particle correlation functions within an Exact Diagonalization (ED) framework.
*   **Role in Project:** Acts as a key computational kernel within an ED solver. It takes an input eigenstate (from an `eigensector_type` object) and an operator specification (creation/destruction, spin, and a mask for sites) and computes the resulting many-body state(s), populating a target `eigensector_type`. The module is equipped to handle different many-body basis representations (Sz-adapted basis or a tensor product of spin-up and spin-down Hilbert spaces) and can leverage MPI for parallel computation in certain scenarios.

## Key Components

*   **Module-Level Variables:**
    *   `INTEGER, POINTER :: IMPiorb(:, :) => NULL()`: (Commented as "not used").
    *   `INTEGER, POINTER :: AIMIMPiorbupdo(:, :) => NULL()`: Stores global orbital indices for impurity sites, mapping (site, spin) to a global index for the up/down product basis.
    *   `INTEGER, POINTER :: AIMIMPiorbsz_updo(:, :) => NULL()`: Stores global orbital indices for impurity sites, mapping (site, spin) to a global index for the Sz-adapted basis. This is used by routines like `apply_Csz_updo`.
    *   `INTEGER, POINTER :: AIMIMPiorbsz(:) => NULL()`: Stores global orbital indices for impurity sites in a flattened way (all up spins then all down spins) for the Sz-adapted basis. Used by `apply_Csz`.
    *   `INTEGER :: Nc`: Number of impurity sites, extracted from the AIM definition.
    *   `LOGICAL, PARAMETER :: force_reset_list = .false.`: A parameter that, if true, forces the deletion of any existing list of eigenstates in the target `eigensector_type` before adding newly computed states.

*   **Public Subroutines:**
    *   **`SUBROUTINE init_apply_C(AIM)`**: Initializes the module by caching orbital mapping information from the provided `AIM_type` object. It populates `Nc`, `AIMIMPiorbupdo`, `AIMIMPiorbsz_updo`, and `AIMIMPiorbsz`. This must be called before other routines.
    *   **`SUBROUTINE Cup_sector(Csec, pm, sector_in)`**: Determines the resulting Hilbert space sector (`Csec`) after applying a spin-up creation (`pm='+'`) or annihilation (`pm='-'`) operator to an initial sector (`sector_in`). It dispatches to either `Csz_sector_updo` or `Cupdo_sector` based on whether `sector_in` is an Sz-adapted (`%sz` is associated) or up/down product basis (`%updo` is associated).
    *   **`SUBROUTINE Cdo_sector(Csec, pm, sector_in)`**: Similar to `Cup_sector`, but for determining the resulting sector after applying a spin-down operator.
    *   **`SUBROUTINE apply_Cup(Ces, pm, MASK, es, esrank)`**: Applies a spin-up operator (creation if `pm='+'`, annihilation if `pm='-'`) to a specific eigenstate `es%lowest%eigen(esrank)`. The `MASK` (logical array) specifies which impurity site(s) the operator acts upon. The resulting state vector(s) are stored in the `Ces` (target `eigensector_type`). It internally calls either `apply_Csz_updo` or `apply_Cupdo`.
    *   **`SUBROUTINE apply_Cdo(Ces, pm, MASK, es, esrank)`**: Similar to `apply_Cup`, but for applying a spin-down operator.
    *   **`SUBROUTINE apply_N_Cup(Ces, pm, MASK, es, esrank)`**: Applies a modified spin-up operator, likely related to Nambu Green's functions or specific correlation functions (e.g., \(N_k c^\dagger_{k,\uparrow}\) or \(N_k c_{k,\uparrow}\)). It currently only supports the Sz-adapted basis (`apply_Csz_N_updo`) and will stop if the basis is up/down product type.
    *   **`SUBROUTINE apply_N_Cdo(Ces, pm, MASK, es, esrank)`**: Similar to `apply_N_Cup`, but for a modified spin-down operator.

*   **Internal Core Subroutines (Private or dispatched by public routines):**
    *   `SUBROUTINE Csz_sector_updo(Csec, pm, spin, sector_in)`: Calculates the target sector for applying an operator of `spin` to an Sz-adapted basis `sector_in`. Note: for `spin=2` (down spin in Nambu context), creation/annihilation roles for `pm` are swapped relative to particle number change.
    *   `SUBROUTINE Cupdo_sector(Csec, pm, spin, sector_in)`: Calculates the target sector for applying an operator of `spin` to an up/down product basis `sector_in`.
    *   `SUBROUTINE apply_Csz_updo(Ces, pm, spin, MASK, es, esrank)`: The workhorse for applying \(c_{site,spin}\) or \(c^\dagger_{site,spin}\) to an eigenstate in the Sz-adapted basis. It iterates over sites specified in `MASK`, retrieves the appropriate global orbital index from `AIMIMPiorbsz_updo`, applies the `create` or `destroy` operation (from `fermion_ket_class`) to each component of the input eigenstate vector, and accumulates the results (including fermion signs) into the target `eigen_out%vec%rc`.
    *   `SUBROUTINE apply_Csz_N_updo(Ces, pm, spin, MASK, es, esrank)`: A specialized version of `apply_Csz_updo` that applies a modified operator. Before applying `create` or `destroy`, it checks the occupation number (`ni__` from `fermion_ket_class`) of related orbitals (e.g., the other spin at the same site for Nambu-like terms).
    *   `SUBROUTINE apply_Cupdo(Ces, pm, spin, MASK, es, esrank)`: The workhorse for applying operators in the up/down product basis. It applies the operator to the relevant spin component of the states. This routine includes logic for MPI parallelism if `USE_TRANSPOSE_TRICK_MPI` is true, using helper subroutines `collect_on_rank0` and `scatter_rank0` to manage distributed eigenvector data.

## Important Variables/Constants

*   **Module-level pointer arrays (orbital mappings):**
    *   `AIMIMPiorbupdo(Nc, 2)`: Stores global orbital indices for impurity sites, mapping `(site_index, spin_index)` to a global orbital index for the up/down product basis representation. `spin_index=1` for up, `spin_index=2` for down.
    *   `AIMIMPiorbsz_updo(Nc, 2)`: Stores global orbital indices for impurity sites, mapping `(site_index, spin_index)` to a global orbital index for the Sz-adapted basis representation.
    *   `AIMIMPiorbsz(Nc*2)`: A flattened version of orbital indices for the Sz-adapted basis (e.g., all spin-up orbitals first, then all spin-down).
    These are initialized by `init_apply_C` and are crucial for ensuring operators act on the correct orbitals in the many-body basis.
*   **`Nc :: INTEGER`**: Number of impurity sites.
*   **`force_reset_list :: LOGICAL, PARAMETER = .false.`**: If set to `.TRUE.` (not default), it would cause `delete_eigenlist(Ces%lowest)` to be called before adding new states generated by operator application. This ensures `Ces` only contains the results of the current operation.

## Usage Examples

```fortran
MODULE example_apply_c_usage
  USE AIM_class, ONLY: AIM_type
  USE apply_C
  USE eigen_sector_class, ONLY: eigensector_type
  USE sector_class, ONLY: sector_type, delete_sector_list ! For cleanup
  USE eigen_class, ONLY: delete_eigenlist ! For cleanup

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_apply_c_module

CONTAINS
  SUBROUTINE test_apply_c_module(aim_model, source_eigenstate, source_rank)
    TYPE(AIM_type), INTENT(IN) :: aim_model
    TYPE(eigensector_type), INTENT(IN) :: source_eigenstate ! Assume this is populated
    INTEGER, INTENT(IN) :: source_rank       ! Rank of the specific eigenvec in source_eigenstate%lowest

    TYPE(eigensector_type) :: result_eigenstate
    LOGICAL, DIMENSION(aim_model%Nc) :: site_operation_mask
    CHARACTER(LEN=1) :: operation_type

    ! Initialize the apply_C module with the AIM orbital information
    CALL init_apply_C(aim_model)

    ! Example: Apply c_dagger_up to the first impurity site (site index 0 in C, 1 in Fortran)
    site_operation_mask = .FALSE.
    site_operation_mask(1) = .TRUE. ! Act on the first impurity site
    operation_type = '+' ! Creation operator

    ! 1. Determine the sector of the resulting state
    CALL Cup_sector(result_eigenstate%sector, operation_type, source_eigenstate%sector)
    ! Note: result_eigenstate%lowest might need initialization if force_reset_list is false
    ! and if the add_eigen routine doesn't handle initially unassociated pointers.
    ! For safety, one might do: CALL new_eigen_list(result_eigenstate%lowest)

    ! 2. Apply the operator
    CALL apply_Cup(Ces=result_eigenstate, pm=operation_type, MASK=site_operation_mask, &
     &             es=source_eigenstate, esrank=source_rank)

    PRINT *, "Applied c_up_dag(site 1). Resulting eigenstate(s) in result_eigenstate."
    ! (result_eigenstate%lowest would now contain the new state vectors)

    ! Cleanup (simplified)
    IF (ASSOCIATED(result_eigenstate%lowest)) THEN
        CALL delete_eigenlist(result_eigenstate%lowest)
    ENDIF
    CALL delete_sector_list(result_eigenstate%sector)


    ! Example: Apply c_down (annihilation) to the second impurity site (if Nc >= 2)
    IF (aim_model%Nc >= 2) THEN
        site_operation_mask = .FALSE.
        site_operation_mask(2) = .TRUE.
        operation_type = '-'

        CALL Cdo_sector(result_eigenstate%sector, operation_type, source_eigenstate%sector)
        CALL apply_Cdo(Ces=result_eigenstate, pm=operation_type, MASK=site_operation_mask, &
         &               es=source_eigenstate, esrank=source_rank)
        PRINT *, "Applied c_down(site 2). Resulting eigenstate(s) in result_eigenstate."
        IF (ASSOCIATED(result_eigenstate%lowest)) THEN
            CALL delete_eigenlist(result_eigenstate%lowest)
        ENDIF
        CALL delete_sector_list(result_eigenstate%sector)
    END IF

  END SUBROUTINE test_apply_c_module
END MODULE example_apply_c_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `log_unit` (logging) and `DP` (double precision kind).
    *   `aim_class`: The `init_apply_C` subroutine requires an `AIM_type` object to set up orbital mappings.
    *   `sector_class`: For `sector_type` definitions and utility functions like `delete_sector`, `npart_func`, `norbs__`, `dimen_func`. This is used to define the Hilbert space sector of the resulting states.
    *   `fermion_hilbert_class`: For `new_fermion_sector` (constructor for Sz-adapted sectors) and `fermion_sector_type`.
    *   `fermion_sector2_class`: For `new_fermion_sector2` (constructor for up/down product basis sectors) and `tabrankupdo` (state ranking in product basis).
    *   `eigen_sector_class`: For the `eigensector_type`, which is the primary data structure for holding input and output eigenstates.
    *   `eigen_class`: For `eigen_type` (representing a single eigenvector and its properties) and its manipulation routines (`new_eigen`, `add_eigen`, `delete_eigen`, `delete_eigenlist`).
    *   `fermion_ket_class`: This is crucial for the low-level operations. Routines like `create` (applies \(c^\dagger\)), `destroy` (applies \(c\)), `new_ket` (creates a ket), and `ni__` (calculates occupation number) are used to act on individual Fock states within an eigenstate's superposition.
    *   `globalvar_ed_solver`: For the `use_transpose_trick_mpi` flag, which controls MPI parallelization strategy in `apply_Cupdo`.
    *   `mpi_mod`: Provides MPI communication primitives like `mpibarrier`, `scatter_it`, `mpigather_on_masternode`, used in `apply_Cupdo` when `use_transpose_trick_mpi` is enabled.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) library, if the parallel features in `apply_Cupdo` are compiled and used.

*   **Interactions with other components:**
    *   This module is a core computational engine within an Exact Diagonalization (ED) solver.
    *   It is initialized by `init_apply_C` using an `AIM_type` object to understand the impurity orbital structure and mapping to the global basis.
    *   It receives `eigensector_type` objects (representing many-body eigenstates) as input for operator application.
    *   The fundamental fermionic operations are delegated to `fermion_ket_class` (e.g., `create`, `destroy` acting on individual kets).
    *   The results (new states after operator application) are stored and returned in an `eigensector_type` object.
    *   For calculations in the up/down product basis (`updo`), the `apply_Cupdo` routine can interact with MPI (via `mpi_mod`) to distribute parts of the calculation across multiple processes, especially if `use_transpose_trick_mpi` is active.
    *   The specialized `apply_N_Cup` and `apply_N_Cdo` routines (and their internal counterparts `apply_Csz_N_updo`) indicate use cases for calculating specific types of correlation functions or Green's functions where operator application is conditional on local occupation numbers (e.g., for Nambu formalism or particular vertex calculations).

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
