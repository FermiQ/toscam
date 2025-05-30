# `src/ed_solver/apply_ns.f90`

## Overview

*   **Purpose:** This module provides a collection of Fortran subroutines designed to apply various local quantum mechanical operators to a given many-body eigenstate. Specifically, it handles the application of number operators (\(N_{site}\)), spin-projection operators (\(S_{z,site}\)), spin raising/lowering operators (\(S^{\pm}_{site}\)), and the Identity operator (\(Id_{site}\)) for specified impurity sites.
*   **Role in Project:** Functions as a computational kernel within an Exact Diagonalization (ED) solver. It takes an input eigenstate (encapsulated in an `eigensector_type` object) and an operator specification (type of operator, site mask) and computes the effect of this operator on the eigenstate. This is crucial for calculating expectation values of these physical observables or for constructing more complex multi-particle correlation functions. The module supports different basis representations (Sz-adapted or Up/Down product basis).

## Key Components

*   **`MODULE apply_NS`**:
    *   **Module-Level Variables:**
        *   `INTEGER, ALLOCATABLE :: IMPiorbupdo(:, :)`: Stores global orbital indices for impurity sites, mapping `(site_index, spin_index)` to a global index. Initialized for up/down product basis context. The spin_index=2 part is noted with a "WARNING!" as it's set equal to spin_index=1 part, which might be specific to how it's used or an oversight if distinct down spin orbitals are needed directly from this variable.
        *   `INTEGER, ALLOCATABLE :: IMPiorbsz(:, :)`: Stores global orbital indices for impurity sites for an Sz-adapted basis, mapping `(site_index, spin_index_in_AIM_IMPiorb)` to a global orbital index.
        *   `INTEGER :: Ns`: Total number of sites in the AIM system (impurity + bath), derived from `AIM%Ns`.
        *   `INTEGER :: Ns2`: Equal to `Ns*2`, likely representing the total number of spin-orbitals in the Sz-adapted basis if `Ns` is the number of spatial sites in that basis.
        *   `LOGICAL, PARAMETER :: force_reset_list = .false.`: If true, forces the deletion of any existing list of eigenstates in the target `eigensector_type` before adding newly computed states.

*   **Public Subroutines:**
    *   **`SUBROUTINE init_apply_NS(AIM)`**: Initializes the module by caching orbital mapping information from the provided `AIM_type` object. It populates `Ns`, `Ns2`, `IMPiorbupdo`, and `IMPiorbsz` using data from `AIM%IMPiorb`. This subroutine must be called before any operator application routines.
    *   **`SUBROUTINE Sz_sector(Csec, pm, sector_in)`**: Determines the resulting sector `Csec` after applying an \(S_z\) operator. Since \(S_z\) is diagonal and does not change particle numbers or total Sz (in an Sz eigenbasis), the output sector `Csec` is a copy of `sector_in`. The `pm` argument is a dummy.
    *   **`SUBROUTINE N_sector(Csec, pm, sector_in)`**: Similar to `Sz_sector`, determines the sector after applying a number operator \(N\). The sector does not change. `pm` is a dummy.
    *   **`SUBROUTINE S_sector(Csec, pm, sector_in)`**: Determines the resulting sector `Csec` after applying \(S^+\) (`pm='+'`) or \(S^-\) (`pm='-'`).
        *   If `sector_in` is Sz-adapted (`%sz` basis): The code `npart_func(sector_in) + 2` for \(S^+\) and `-2` for \(S^-\) suggests these are not standard S+ or S- operators (which conserve total particle number). This might implement pair creation (\(c^\dagger_\uparrow c^\dagger_\downarrow\)) for \(S^+\) and pair annihilation (\(c_\downarrow c_\uparrow\)) for \(S^-\).
        *   If `sector_in` is an up/down product basis (`%updo` basis): \(S^+\) changes `(Nup, Ndo)` to `(Nup+1, Ndo-1)`. \(S^-\) changes `(Nup, Ndo)` to `(Nup-1, Ndo+1)`. This is the standard definition for \(S^\pm\) operators.
    *   **`SUBROUTINE apply_Sz(Ces, pm, MASK, es, esrank)`**: Applies the local spin-projection operator \(S_{z,site}\) to the eigenstate `es%lowest%eigen(esrank)` for sites specified in `MASK`.
        *   In Sz-basis: Calculates \(S_{z,site} = 0.5 \times (n_{site,\uparrow} - (1 - n_{site,\downarrow}))\). This Nambu-like definition suggests \( (1 - n_{site,\downarrow}) \) represents the number of holes with spin down.
        *   In up/down product basis: Calculates \(S_{z,site} = 0.5 \times (n_{site,up} - n_{site,down})\).
        Since \(S_z\) is diagonal, this routine multiplies the components of the input eigenvector by their respective \(S_z\) eigenvalues.
    *   **`SUBROUTINE apply_S(Ces, pm, MASK, es, esrank)`**: Applies local spin raising \(S^+_{site}\) (if `pm='+'`) or lowering \(S^-_{site}\) (if `pm='-'`) operators.
        *   In Sz-basis: `pm='+'` corresponds to `create_pair(IMPiorbsz(site,1), IMPiorbsz(site,2), ket_in)`, which is \(c^\dagger_{site,\uparrow} c^\dagger_{site,\downarrow}\). `pm='-'` corresponds to `destroy_pair(IMPiorbsz(site,2), IMPiorbsz(site,1), ket_in)`, which is \(c_{site,\downarrow} c_{site,\uparrow}\). These are pair creation/annihilation operators, not standard S+/S-.
        *   In up/down product basis: `pm='+'` applies \(c^\dagger_{site,up} c_{site,down}\). `pm='-'` applies \(-c^\dagger_{site,down} c_{site,up}\) (standard S+ and S- definitions, with the conventional minus sign for S-).
    *   **`SUBROUTINE apply_N(Ces, pm, MASK, es, esrank)`**: Applies the local number operator \(N_{site}\) for sites specified in `MASK`.
        *   In Sz-basis: Calculates \(N_{site} = n_{site,\uparrow} + (1 - n_{site,\downarrow})\) (Nambu particle number).
        *   In up/down product basis: Calculates \(N_{site} = n_{site,up} + n_{site,down}\).
        This is a diagonal operator.
    *   **`SUBROUTINE apply_Id(Ces, pm, MASK, es, esrank)`**: Applies the Identity operator for the sites specified in `MASK`. This effectively copies the input eigenstate `es%lowest%eigen(esrank)` to the output `Ces`, associating it with the specified `site` rank in `Ces%lowest`. The `pm` argument is a dummy.

## Important Variables/Constants

*   **Module-level pointer arrays for orbital mappings (initialized by `init_apply_NS`):**
    *   `IMPiorbupdo(AIM%Nc, 2)`: Stores global orbital indices for impurity sites in the context of an up/down product basis. `IMPiorbupdo(site,1)` for up-spin orbital, `IMPiorbupdo(site,2)` is set to the same as spin-up (marked with "WARNING!").
    *   `IMPiorbsz(AIM%Nc, 2)`: Stores global orbital indices for impurity sites for an Sz-adapted basis. `IMPiorbsz(site,1)` for the up-spin component and `IMPiorbsz(site,2)` for the down-spin component mapping from the original AIM definition.
*   **`Ns :: INTEGER`**: Total number of spatial sites in the AIM system (impurity + bath).
*   **`Ns2 :: INTEGER`**: `Ns * 2`, likely representing the total number of spin-orbitals if `Ns` is the number of sites in the Sz-adapted basis.
*   **`force_reset_list :: LOGICAL, PARAMETER = .false.`**: If true, the list of eigenstates in the target `Ces` is cleared before adding new states resulting from the operator application.

## Usage Examples

```fortran
MODULE example_apply_ns_usage
  USE AIM_class, ONLY: AIM_type
  USE apply_NS
  USE eigen_sector_class, ONLY: eigensector_type
  USE sector_class, ONLY: sector_type, delete_sector_list
  USE eigen_class, ONLY: delete_eigenlist

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_apply_ns_module

CONTAINS
  SUBROUTINE test_apply_ns_module(aim_model, source_eigenstate, source_rank)
    TYPE(AIM_type), INTENT(IN) :: aim_model
    TYPE(eigensector_type), INTENT(IN) :: source_eigenstate ! Assume populated
    INTEGER, INTENT(IN) :: source_rank       ! Rank of specific eigenvec in source_eigenstate

    TYPE(eigensector_type) :: result_eigenstate
    LOGICAL, DIMENSION(aim_model%Nc) :: site_mask_example

    ! Initialize the apply_NS module
    CALL init_apply_NS(aim_model)

    ! Example: Apply Sz operator to the first impurity site
    site_mask_example = .FALSE.
    site_mask_example(1) = .TRUE.

    CALL Sz_sector(result_eigenstate%sector, pm=' ', sector_in=source_eigenstate%sector)
    CALL apply_Sz(Ces=result_eigenstate, pm=' ', MASK=site_mask_example, &
     &            es=source_eigenstate, esrank=source_rank)
    ! result_eigenstate now holds Sz_site1 * |source_eigenstate>
    ! (Operator is diagonal, so only coefficients change, or become zero)
    PRINT *, "Applied Sz to site 1."
    IF (ASSOCIATED(result_eigenstate%lowest)) CALL delete_eigenlist(result_eigenstate%lowest)
    CALL delete_sector_list(result_eigenstate%sector)


    ! Example: Apply N operator to the first impurity site
    CALL N_sector(result_eigenstate%sector, pm=' ', sector_in=source_eigenstate%sector)
    CALL apply_N(Ces=result_eigenstate, pm=' ', MASK=site_mask_example, &
     &           es=source_eigenstate, esrank=source_rank)
    PRINT *, "Applied N to site 1."
    IF (ASSOCIATED(result_eigenstate%lowest)) CALL delete_eigenlist(result_eigenstate%lowest)
    CALL delete_sector_list(result_eigenstate%sector)


    ! Example: Apply S+ operator (pair creation in Sz basis, c_up_dag*c_down in updo basis)
    ! to the first impurity site
    CALL S_sector(result_eigenstate%sector, pm='+', sector_in=source_eigenstate%sector)
    CALL apply_S(Ces=result_eigenstate, pm='+', MASK=site_mask_example, &
     &           es=source_eigenstate, esrank=source_rank)
    PRINT *, "Applied S+ to site 1."
    IF (ASSOCIATED(result_eigenstate%lowest)) CALL delete_eigenlist(result_eigenstate%lowest)
    CALL delete_sector_list(result_eigenstate%sector)

  END SUBROUTINE test_apply_ns_module
END MODULE example_apply_ns_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP` (double precision kind).
    *   `aim_class`: The `init_apply_NS` subroutine requires an `AIM_type` object to set up orbital mappings.
    *   `sector_class`: For `sector_type` and its utility functions (`delete_sector`, `new_sector`, `npart_func`, `dimen_func`).
    *   `fermion_hilbert_class`: For `new_fermion_sector` (used in `S_sector` for Sz-basis).
    *   `fermion_sector2_class`: For `new_fermion_sector2` (used in `S_sector` for up/down basis) and `rankupdo` (used in `apply_Sz`, `apply_S`, `apply_N` for up/down basis state indexing).
    *   `eigen_sector_class`: For `eigensector_type`, which encapsulates eigenstates.
    *   `eigen_class`: For `eigen_type` (representing a single eigenstate vector) and its manipulation routines (`add_eigen`, `delete_eigen`, `delete_eigenlist`, `new_eigen`).
    *   `fermion_ket_class`: For `fermion_ket_type`, `new_ket` (to create kets from integer representations), `ni__` (to compute occupation number of an orbital in a ket), `create_pair`, `destroy_pair`, `create`, `destroy`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules. No explicit third-party external libraries are directly used.

*   **Interactions with other components:**
    *   This module is a utility within an ED solver, used for calculating the effect of various local physical operators.
    *   It is initialized using `AIM_type` data to correctly map local site and spin operations to the global orbital indices of the AIM.
    *   It operates on `eigensector_type` objects, which are expected to be populated by diagonalization routines.
    *   The core operations of applying operators to individual Fock states (kets) are delegated to routines from `fermion_ket_class`.
    *   The definitions of \(S_z\), \(N\), and particularly \(S^\pm\) in the Sz-adapted basis (using `create_pair`/`destroy_pair` and Nambu-like number definitions) suggest specialized use cases, possibly related to superconductivity or models where such pair operators are relevant. In the standard up/down product basis, \(S_z\), \(N\), \(S^\pm\) follow more conventional definitions.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
