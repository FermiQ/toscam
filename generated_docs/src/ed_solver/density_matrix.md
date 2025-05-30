# `src/ed_solver/density_matrix.f90`

## Overview

*   **Purpose:** This module is dedicated to computing the reduced density matrix for the impurity subsystem of an Anderson Impurity Model (AIM). It also provides routines to analyze this density matrix, such as calculating the von Neumann entanglement entropy and identifying the dominant impurity eigenstates and their occupation probabilities.
*   **Role in Project:** Serves as an important post-processing and analysis tool within an Exact Diagonalization (ED) solver framework. After the full AIM Hamiltonian is diagonalized, this module allows for focusing on the properties of the impurity degrees of freedom by tracing out the bath. The resulting reduced density matrix provides complete information about the quantum state of the impurity, including entanglement with the bath.

## Key Components

*   **`MODULE density_matrix`**:
    The main module encapsulating all related subroutines.

*   **Public Subroutines:**
    *   **`SUBROUTINE compute_density_matrix(dmat, AIM, beta, GS)`**:
        Calculates the reduced density matrix (`dmat`) for the impurity sites.
        *   **Arguments:**
            *   `dmat (TYPE(rcmatrix_type), INTENT(INOUT))`: The output impurity reduced density matrix. `rcmatrix_type` likely stores a real or complex matrix of size `nIMPstates x nIMPstates`.
            *   `AIM (TYPE(AIM_type), INTENT(IN))`: The definition of the full Anderson Impurity Model, providing context about impurity and bath Hilbert spaces.
            *   `beta (REAL(DP), INTENT(IN))`: The inverse temperature \(1/(k_B T)\). Used for Boltzmann weighting of eigenstates.
            *   `GS (TYPE(eigensectorlist_type), INTENT(IN))`: The list of low-lying many-body eigenstates (and eigenvalues) of the full AIM, obtained from an ED calculation.
        *   **Functionality:**
            1.  Initializes the output density matrix `dmat%rc` to zero.
            2.  Determines the total number of impurity states (`nIMPstates`) and bath states (`nBATHstates`).
            3.  Calculates the partition function `Zpart = \sum_k e^{-\beta (E_k - E_0)}` and the ground state energy `E0` from `GS`.
            4.  The computation of `dmat` elements \(\rho_{ij} = \langle i | \hat{\rho}_{imp} | j \rangle\) involves a sum over all eigenstates \(| \Psi_k \rangle\) of the full system:
                \[ \rho_{ij} = \frac{1}{Z_{part}} \sum_k e^{-\beta (E_k - E_0)} \sum_{bath\_state} \langle i, bath\_state | \Psi_k \rangle \langle \Psi_k | j, bath\_state \rangle \]
            5.  This is implemented by iterating through each eigenstate `eigen` in each sector `es` of `GS`.
            6.  For each `eigen`, its Boltzmann factor `boltz = exp(-beta*(eigen%val - E0))` is computed.
            7.  The code then loops over pairs of impurity many-body basis states (`IMPrank1`, `IMPrank2`). This part is parallelized using OpenMP (`!$OMP PARALLEL DO`) and can be distributed across MPI processes.
            8.  Inside this loop, it performs a trace over all bath basis states (`BATHrank`):
                *   For a given `IMPrank1` (impurity state \(|i\rangle\)), `IMPrank2` (impurity state \(|j\rangle\)), and `BATHrank` (bath state \(|b\rangle\)), it constructs the corresponding full AIM state indices `AIMstate1` (for \(|i,b\rangle\)) and `AIMstate2` (for \(|j,b\rangle\)) using `IMPBATH2AIMstate` (from `aim_class.f90`).
                *   It retrieves the amplitudes of these full AIM states in the current eigenstate: `coeff1 = <i,b | eigen>` and `coeff2 = <j,b | eigen>`.
                *   The contribution `coeff1 * CONJG(coeff2) * boltz` is added to `dmat%rc(IMPrank1, IMPrank2)`.
                *   Hermiticity `dmat%rc(IMPrank2, IMPrank1) = dmat%rc(IMPrank2, IMPrank1) + CONJG(coeff1) * coeff2 * boltz` is ensured if `IMPrank1 /= IMPrank2`.
            9.  If running with MPI (`size2 > 1`), an internal subroutine `mpi_collect_` is called to sum up the contributions to `dmat` from all processes.
            10. Finally, the computed `dmat` is normalized by the partition function `Zpart`.

    *   **`SUBROUTINE analyze_density_matrix(dmat, IMPiorb, NAMBU)`**:
        Analyzes the computed impurity reduced density matrix `dmat`.
        *   **Arguments:**
            *   `dmat (REAL(DP) or COMPLEX(DP), INTENT(IN))`: The impurity reduced density matrix (passed as `dmat(:,:)` which is `dmat%rc%mat` if `dmat` is `rcmatrix_type`).
            *   `IMPiorb (INTEGER, INTENT(IN))`: The impurity orbital mapping array (from `AIM_type%impurity%iorb`), used by `cket_from_state` to interpret the impurity basis states.
            *   `NAMBU (LOGICAL, INTENT(IN))`: A flag indicating whether the impurity states are in a Nambu (superconducting) basis, affecting how `cket_from_state` formats them.
        *   **Functionality:**
            1.  Diagonalizes the density matrix `dmat` using `diagonalize` (from `matrix` module) to obtain its eigenvalues (`VALP`, which are occupation probabilities) and eigenvectors (`VECP`).
            2.  Prints the `ed_num_eigenstates_print` largest eigenvalues and their corresponding eigenvectors. Eigenvectors are displayed as linear combinations of impurity basis states, with the basis states formatted by `cket_from_state` (from `readable_vec_class`).
            3.  Calculates the von Neumann entanglement entropy: \(S = - \sum_i \lambda_i \log(\lambda_i)\), where \(\lambda_i\) are the eigenvalues (`VALP`) of `dmat`. This quantifies the entanglement between the impurity and the bath.
            4.  Computes and prints the trace of the density matrix (Tr(\(\rho\))), which should be close to 1 for a correctly normalized density matrix.

*   **Internal Subroutine `mpi_collect_` (within `compute_density_matrix`)**:
    This subroutine is responsible for gathering the distributed parts of the density matrix that were computed by different MPI processes. It uses `MPI_ALLGATHERV` to collect the lower (or upper, due to symmetry) triangular parts of `dmat` from each process and assembles the full matrix on all processes, ensuring Hermiticity.

## Important Variables/Constants

*   **`dmat :: TYPE(rcmatrix_type)`**: The primary data structure holding the impurity reduced density matrix. It's an output of `compute_density_matrix` and an input to `analyze_density_matrix`.
*   **`AIM :: TYPE(AIM_type)`**: Provides the definitions of the impurity and bath Hilbert spaces and their mapping, essential for `compute_density_matrix`.
*   **`GS :: TYPE(eigensectorlist_type)`**: Contains the many-body eigenstates and eigenvalues of the full AIM system, which are the starting point for calculating the density matrix.
*   **`beta :: REAL(DP)`**: The inverse temperature, used for Boltzmann weighting of the eigenstates when constructing the thermal density matrix.
*   **`VALP(:) :: REAL(DP)` (in `analyze_density_matrix`)**: Array storing the eigenvalues (occupation probabilities) of the reduced density matrix.
*   **`VECP(:,:)` (in `analyze_density_matrix`)**: Matrix storing the eigenvectors of the reduced density matrix. Each column is an eigenvector.
*   **`ENTROPY :: REAL(DP)` (in `analyze_density_matrix`)**: The calculated von Neumann entanglement entropy of the impurity subsystem.

## Usage Examples

```fortran
MODULE example_density_matrix_usage
  USE density_matrix
  USE aim_class, ONLY: AIM_type
  USE impurity_class, ONLY: impurity_type ! For AIM%impurity%iorb
  USE eigen_sector_class, ONLY: eigensectorlist_type
  USE rcmatrix_class, ONLY: rcmatrix_type, new_rcmatrix, delete_rcmatrix
  USE genvar, ONLY: DP
  ! ... other necessary modules for AIM and ED solution setup ...

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: test_density_matrix_computation

CONTAINS
  SUBROUTINE test_density_matrix_computation
    TYPE(AIM_type) :: aim_model_instance
    TYPE(eigensectorlist_type) :: ed_solution_gs
    REAL(DP) :: inverse_temp
    TYPE(rcmatrix_type) :: impurity_rho

    ! --- Assume aim_model_instance and ed_solution_gs are properly initialized ---
    ! (e.g., aim_model_instance%impurity%nstates is set)
    ! inverse_temp = 50.0 ! Example value
    ! --- End of initialization ---

    ! Allocate and initialize the density matrix structure
    CALL new_rcmatrix(impurity_rho, "ImpurityRho", &
     &                aim_model_instance%impurity%nstates, aim_model_instance%impurity%nstates)

    PRINT *, "Computing impurity reduced density matrix..."
    CALL compute_density_matrix(impurity_rho, aim_model_instance, inverse_temp, ed_solution_gs)
    PRINT *, "Density matrix computation finished."

    ! Analyze the computed density matrix
    PRINT *, "Analyzing density matrix..."
    CALL analyze_density_matrix(impurity_rho%rc%mat, aim_model_instance%impurity%iorb, &
     &                          NAMBU=aim_model_instance%bath%SUPER) ! Use SUPER flag from bath
                                                                     ! to indicate Nambu basis

    ! Clean up
    CALL delete_rcmatrix(impurity_rho)
    ! (Cleanup of aim_model_instance, ed_solution_gs would happen here)

  END SUBROUTINE test_density_matrix_computation
END MODULE example_density_matrix_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `eigen_sector_class`: For `eigensector_type`, `eigensectorlist_type` (input for eigenstates), and utility functions `gsenergy`, `partition`.
    *   `genvar`: For `DP` (double precision kind) and MPI-related variables (`iproc`, `no_mpi`, `nproc`, `size2`, `ierr`). Also `log_unit`.
    *   `linalg`: For `conj` (complex conjugate) and `dexpc` (complex exponential `exp(-beta*E)`).
    *   `common_def`: For messaging (`dump_message`), timing (`reset_timer`, `timer_fortran`), and string utilities (`c2s`, `i2c`).
    *   `sector_class`: For `is_in_sector` (to check if a constructed AIM state belongs to the current eigenstate's sector) and `rank_func` (to get the index of an AIM state within a sector).
    *   `aim_class`: For `AIM_type` (input model definition) and `IMPBATH2AIMstate` (to map impurity/bath state pairs to full AIM state indices).
    *   `mpi_mod`: For `split` (to distribute work among MPI processes).
    *   `rcmatrix_class`: For `rcmatrix_type`, the data structure used to store the density matrix.
    *   `eigen_class`: For `eigen_type`, representing a single eigenstate vector.
    *   `timer_mod`: For `start_timer`, `stop_timer`.
    *   `matrix` (used in `analyze_density_matrix`): For the `diagonalize` routine and `write_array`.
    *   `readable_vec_class` (used in `analyze_density_matrix`): For `cket_from_state`, used to print eigenvectors in a human-readable Fock state basis.
    *   `sorting` (used in `analyze_density_matrix`): For `qsort_array` to sort eigenvector components by magnitude.

*   **External Libraries:**
    *   **MPI (Message Passing Interface):** The `compute_density_matrix` subroutine uses MPI for parallel computation, specifically `MPI_ALLGATHERV` in the internal `mpi_collect_` routine to sum contributions from different processes. The standard Fortran `mpi` module is used.

*   **Interactions with other components:**
    *   This module is a post-processing tool that operates on the results of an ED calculation.
    *   It takes the Anderson Impurity Model definition (`AIM_type`) and the list of its computed eigenstates and eigenvalues (`eigensectorlist_type`) as primary inputs.
    *   The `compute_density_matrix` routine performs a trace over the bath degrees of freedom. This involves:
        *   Iterating through the eigenstates provided in `GS`.
        *   For each eigenstate, iterating through all combinations of impurity basis states and bath basis states.
        *   Using `IMPBATH2AIMstate` (from `aim_class`) to find the representation of these product states in the full AIM basis used for the eigenstate vector.
        *   Summing up the appropriate outer products of eigenvector coefficients, weighted by Boltzmann factors.
    *   The resulting impurity reduced density matrix (`dmat`) can then be used for various purposes, such as calculating expectation values of impurity-only operators or, as done in `analyze_density_matrix`, calculating entanglement measures like the von Neumann entropy.
    *   `analyze_density_matrix` further processes `dmat` by diagonalizing it to find impurity occupation probabilities and the impurity's entanglement with the bath.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
