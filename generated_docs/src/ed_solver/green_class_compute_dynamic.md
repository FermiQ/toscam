# `src/ed_solver/green_class_compute_dynamic.f90`

## Overview

*   **Purpose:** This module provides subroutines to compute the dynamic (frequency-dependent) part of Green's functions or other two-time correlation functions. It supports various numerical methods, primarily Lanczos iteration for an efficient calculation of spectral functions, and direct summation over all eigenstates when a full exact diagonalization (Full ED) solution is available.
*   **Role in Project:** Serves as a computational engine for calculating dynamic response functions within an Exact Diagonalization (ED) solver framework. It takes an initial state (typically an operator applied to the ground state, \(O^\dagger|\Psi_0\rangle\)), its energy, a frequency mesh, and other parameters, then computes the corresponding Green's function component \(G(\omega) = \langle \Psi_0 | O ( \omega - (H-E_0) )^{-1} O^\dagger | \Psi_0 \rangle\) (for particle part) or \(G(\omega) = \langle \Psi_0 | O^\dagger ( \omega + (H-E_0) )^{-1} O | \Psi_0 \rangle\) (for hole part). The module can dispatch to different algorithms based on global settings.

## Key Components

*   **`MODULE green_class_compute_dynamic`**:
    The main module containing the public dispatcher routine and the underlying computational subroutines.

*   **Public Subroutines:**
    *   **`SUBROUTINE compute_dynamic(iph, dyn, freq, Opm, stat, title, normvec, iisector, GS, keldysh_level)`**:
        This is the primary public interface and acts as a dispatcher to the appropriate computational routine based on the global variable `which_lanczos` (from `globalvar_ed_solver.f90`).
        *   **Arguments:**
            *   `iph (INTEGER, INTENT(IN))`: Specifies the part of the spectrum: 1 for particle part (poles at \(\omega > 0\) relative to chemical potential, involves \((z - (H-E_0))^{-1}\)), 2 for hole part (poles at \(\omega < 0\), involves \((z + (H-E_0))^{-1}\)).
            *   `dyn (COMPLEX(DP), INTENT(INOUT))`: A 1D array to store the computed dynamic function values at each frequency point in `freq`.
            *   `freq (TYPE(freq_type), INTENT(IN))`: The frequency mesh (real or Matsubara) on which the function is computed.
            *   `Opm (TYPE(eigen_type), INTENT(IN))`: Represents the initial state for the Green's function calculation. `Opm%vec%rc` is the state vector (e.g., \(O^\dagger|\Psi_0\rangle\)), and `Opm%val` is its reference energy (e.g., \(E_0\)).
            *   `stat (CHARACTER(LEN=*), INTENT(IN))`: Particle statistics, typically `FERMIONIC` or `BOSONIC`, which determines sign factors.
            *   `title (CHARACTER(LEN=*), INTENT(IN))`: A title for logging and output related to this specific calculation.
            *   `normvec (REAL(DP), INTENT(IN))`: The precomputed norm of `Opm%vec%rc`. If close to zero, computation might be skipped.
            *   `iisector (INTEGER, OPTIONAL, INTENT(IN))`: Index of the relevant sector in `GS`. Used only by `compute_dynamic_full_ed`.
            *   `GS (TYPE(eigensectorlist_type), OPTIONAL, INTENT(IN))`: The list of all computed eigenstates and eigenvalues from a full diagonalization. Used only by `compute_dynamic_full_ed`.
            *   `keldysh_level (INTEGER, OPTIONAL, INTENT(IN))`: If present, indicates a Keldysh calculation, currently dispatched to `compute_dynamic_cpu`.
        *   **Dispatch Logic:**
            *   If `keldysh_level` is present, always calls `compute_dynamic_cpu`.
            *   Based on `which_lanczos`:
                *   `'NORMAL'` or `'ARPACK'`: Calls `compute_dynamic_cpu` (standard Lanczos method).
                *   `'GPU'`: (Currently commented out) Would call a GPU-accelerated Lanczos version.
                *   `'FULL_ED'`: If `FLAG_FULL_ED_GREEN` (global flag) is `.FALSE.`, it falls back to `compute_dynamic_cpu`. If `.TRUE.`, it calls `compute_dynamic_full_ed`.
            *   Stops if `which_lanczos` is not recognized.

*   **Internal Computational Subroutines:**
    *   **`SUBROUTINE compute_dynamic_full_ed(iph, dyn, freq, Opm, stat, title, normvec, iisector, GS)`**:
        Computes the Green's function directly using its Lehmann representation by summing over all available eigenstates.
        *   Formula for \(G(\omega)\) components: \(\sum_k \frac{|\langle k | \text{OpState} \rangle|^2}{\omega - \text{ph_sign} \cdot (E_k - E_{ref})} \cdot (-\text{fermion_sign})\), where `OpState` is `Opm%vec%rc`, \(E_{ref}\) is `Opm%val`. The sum is over all states \(|k\rangle\) in `GS%es(iisector)%lowest`.
        *   `norm_ = |\langle k | OpState \rangle|^2` is computed using `MPI_DOT_PRODUCT` or `DOT_PRODUCT` depending on MPI usage.
        *   The loop over eigenstates `iter` is parallelized across MPI processes (`DO iter = rank + 1, ..., size2`). Results are summed using `mpisum`.
    *   **`SUBROUTINE compute_dynamic_cpu(iph, dyn, freq, Opm, stat, title, normvec, keldysh_level)`**:
        Computes the Green's function using a Lanczos iteration method.
        *   Initializes a Lanczos recursion with the starting vector `initvec%rc = Opm%vec%rc / sqrt(norm_Opmvec2)`.
        *   Iteratively calls `one_step_Lanczos_fast` (from `lanczos_fast.f90`) to generate the tridiagonal Lanczos matrix \(T\) (coefficients \(\alpha_i, \beta_i\)), storing it in `Lanczos_matrix`.
        *   The Green's function is then obtained from the continued fraction form: \(G(z) = \frac{\langle \text{initvec} | \text{initvec} \rangle}{z - \alpha_1 - \frac{\beta_2^2}{z - \alpha_2 - \frac{\beta_3^2}{\dots}}}\). This is calculated as `norm_Opmvec2 * (-fermion_sign) * invert_zmtridiag(freq%vec(iw) - ph_sign*Opm%val, tri)`, where `tri` is the (shifted) Lanczos matrix.
        *   If `keldysh_level` is present, it performs a more complex calculation:
            *   The Lanczos matrix (up to `Niter_` iterations) is diagonalized.
            *   The initial vector and Lanczos basis vectors are transformed into this eigenbasis.
            *   A diagonal matrix `DD(iter,iter) = exp(-VALP(iter)*keldysh_delta)` is constructed.
            *   The final state for Keldysh evolution is computed by transforming back, effectively calculating \(e^{-H \delta t} |\text{initvec}\rangle\). The result overwrites `Opm%vec%rc`. (The actual Green's function for Keldysh is not directly computed here, rather the evolved state).
    *   **Internal Functions `internal_scalprod_` and `internal_scalprod` (within `compute_dynamic_cpu`)**: Custom complex dot product functions.

## Important Variables/Constants

*   **Input `Opm :: TYPE(eigen_type)`**: This is a key input. `Opm%vec%rc` provides the initial vector for the Lanczos recursion or the state \(O^\dagger|\Psi_0\rangle\) for the Lehmann sum. `Opm%val` provides the reference energy \(E_0\).
*   **Input `freq :: TYPE(freq_type)`**: Defines the frequency mesh (real or Matsubara points) at which the dynamic function `dyn` will be computed.
*   **Output `dyn(:) :: COMPLEX(DP)`**: The array where the computed Green's function values at each frequency point are stored.
*   **Global `which_lanczos :: CHARACTER(20)` (from `globalvar_ed_solver.f90`)**: Controls which computational method (`cpu` Lanczos, `FULL_ED`, or formerly `GPU` Lanczos) is dispatched.
*   **Global `FLAG_FULL_ED_GREEN :: LOGICAL` (from `globalvar_ed_solver.f90`)**: If `which_lanczos == 'FULL_ED'`, this flag determines if the full ED summation is used or if it falls back to `compute_dynamic_cpu`.
*   `iph :: INTEGER`: Determines particle (1) or hole (2) part of the spectrum, influencing energy denominators.
*   `stat :: CHARACTER(LEN=*)`: `FERMIONIC` or `BOSONIC`, determining overall sign factors.

## Usage Examples

```fortran
MODULE example_compute_dynamic_usage
  USE green_class_compute_dynamic
  USE eigen_class, ONLY: eigen_type
  USE frequency_class, ONLY: freq_type
  USE eigen_sector_class, ONLY: eigensectorlist_type ! For FULL_ED case
  USE genvar, ONLY: DP, FERMIONIC
  USE globalvar_ed_solver, ONLY: which_lanczos, FLAG_FULL_ED_GREEN ! To control dispatch

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: calculate_example_greens_function

CONTAINS
  SUBROUTINE calculate_example_greens_function
    TYPE(eigen_type) :: O_dagger_psi0_state
    TYPE(freq_type) :: frequency_mesh
    COMPLEX(DP), ALLOCATABLE :: G_omega(:)
    REAL(DP) :: norm_O_dagger_psi0
    ! For FULL_ED example:
    TYPE(eigensectorlist_type) :: all_system_eigenstates
    INTEGER :: relevant_sector_idx

    ! --- Assume O_dagger_psi0_state, frequency_mesh are initialized ---
    ! --- Assume norm_O_dagger_psi0 is pre-calculated ---
    ! --- For FULL_ED, all_system_eigenstates and relevant_sector_idx also initialized ---
    ! Example:
    ! ALLOCATE(O_dagger_psi0_state%vec%rc(100)); O_dagger_psi0_state%vec%rc = ...
    ! O_dagger_psi0_state%val = -10.0_DP ! Ground state energy E0
    ! norm_O_dagger_psi0 = 1.0_DP
    ! CALL new_freq(frequency_mesh, Nw=200, wmin=-5.0_DP, wmax=5.0_DP, width=0.1_DP, sense="RETARDED")
    ! ALLOCATE(G_omega(frequency_mesh%Nw))

    ! Set the desired computation method via global variable
    which_lanczos = 'NORMAL' ! or 'FULL_ED'
    FLAG_FULL_ED_GREEN = .TRUE. ! if which_lanczos == 'FULL_ED'

    CALL compute_dynamic(iph=1, dyn=G_omega, freq=frequency_mesh, Opm=O_dagger_psi0_state, &
     &                   stat=FERMIONIC, title="Example_G_up_particle", &
     &                   normvec=norm_O_dagger_psi0, &
     &                   iisector=relevant_sector_idx, GS=all_system_eigenstates)
                         ! iisector and GS only strictly needed if which_lanczos='FULL_ED'
                         ! and FLAG_FULL_ED_GREEN=.TRUE.

    PRINT *, "Computed Green's function G_omega(1) = ", G_omega(1)

    ! --- Cleanup ---
    ! DEALLOCATE(O_dagger_psi0_state%vec%rc)
    ! CALL delete_frequency(frequency_mesh)
    ! DEALLOCATE(G_omega)
    ! CALL delete_eigensectorlist(all_system_eigenstates)

  END SUBROUTINE calculate_example_greens_function
END MODULE example_compute_dynamic_usage
```

## Dependencies and Interactions

*   **Internal Dependencies (within the same project/library):**
    *   `genvar`: For `DP`, `BOSONIC`, `FERMIONIC`, MPI variables (`rank`, `size2`), logging (`log_unit`).
    *   `eigen_class`: For `eigen_type`.
    *   `eigen_sector_class`: For `eigensectorlist_type` (used by `compute_dynamic_full_ed`).
    *   `frequency_class`: For `freq_type`.
    *   `globalvar_ed_solver`: For various global control flags like `FLAG_FULL_ED_GREEN`, `which_lanczos`, `cutoff_dynamic`, `nitergreenmax`, `USE_TRANSPOSE_TRICK_MPI`, `keldysh_delta`, `keldysh_ortho`.
    *   `common_def`: For `utils_assert`.
    *   `h_class`: For `dimen_H()` (Hamiltonian dimension) and `sector_H` (current sector, for `USE_TRANSPOSE_TRICK_MPI`).
    *   `lanczos_fast`: The `one_step_Lanczos_fast` routine is the core of the Lanczos iteration in `compute_dynamic_cpu`.
    *   `linalg`: For `norme` (vector norm), `MPLX` (complex conversion).
    *   `mpi_mod`: For `mpi_dot_product`, `mpibarrier`, `mpisum`.
    *   `rcvector_class`: For `rcvector_type` and its utilities (`delete_rcvector`, `new_rcvector`, `norm_rcvector`).
    *   `tridiag_class`: For `tridiag_type` (to store Lanczos matrix) and its methods (`delete_tridiag`, `diagonalize_tridiag`, `invert_zmtridiag`, `new_tridiag`, `submatrix_tridiag`).
    *   `stringmanip`: For `tostring`.
    *   `timer_mod`: For `start_timer`, `stop_timer`.

*   **External Libraries:**
    *   Standard Fortran intrinsic modules.
    *   MPI (Message Passing Interface) if parallel features are enabled and used (e.g., in `compute_dynamic_full_ed` or potentially within Lanczos if `USE_TRANSPOSE_TRICK_MPI` implies distributed `Hmult`).

*   **Interactions with other components:**
    *   This module is a core computational engine, typically called by routines in `correlations.f90` (specifically, by procedures likely located in the included `correlations_dyn_stat.h` file).
    *   It takes an initial state vector `Opm%vec%rc` (which is \(O^\dagger |\Psi_0\rangle\), where \(O^\dagger\) is a creation-like operator and \(|\Psi_0\rangle\) is often the ground state) and its reference energy `Opm%val` (often \(E_0\)). This initial state is prepared by modules like `apply_c.f90` or `apply_ns.f90`.
    *   The Hamiltonian's action (\(H \cdot \text{vector}\)) is provided by `h_class.f90` through `one_step_Lanczos_fast` which calls `Hmult__`.
    *   If the "FULL_ED" method is chosen, it directly uses the list of all eigenstates (`GS` of `eigensectorlist_type`) obtained from a full diagonalization routine (like `ED_diago` in `ed_arpack.f90`).
    *   The computed dynamic array `dyn` is then used by the calling routines in `correlations.f90` to populate the `fctn` component of `correl_type` or `green_type` objects.

---
_This document is autogenerated. Manual edits may be required for accuracy and completeness._
