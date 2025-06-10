# KSP (Krylov Subspace Methods)

## Overview

The KSP (Krylov Subspace Methods) component in PETSc provides a robust and flexible framework for solving linear systems of equations, typically of the form $Ax = b$. It encompasses a wide variety of iterative methods that construct approximate solutions within a Krylov subspace. KSP solvers are fundamental to many PETSc applications, serving as the linear algebra backbone for more complex solvers like SNES (for nonlinear systems) and TS (for time-dependent problems). The `src/ksp/ksp/interface/itfunc.c` file, along with its primary header `petscksp.h`, defines the main user interface for creating, configuring, and applying these Krylov solvers.

## Key Components (Functions)

The KSP interface includes functions for managing the lifecycle and settings of linear solvers:

*   **`KSPCreate(MPI_Comm comm, KSP *ksp)`**: Creates a KSP context.
    *   Input: `comm` - MPI communicator.
    *   Output: `ksp` - The newly created KSP object.
*   **`KSPSetType(KSP ksp, KSPType type)`**: Sets the specific Krylov method to be used (e.g., `KSPGMRES`, `KSPCG`).
    *   Input: `ksp` - The KSP context.
    *   Input: `type` - A string representing the desired KSP type.
*   **`KSPSetOperators(KSP ksp, Mat Amat, Mat Pmat)`**: Sets the operator matrix (A) that defines the linear system and the matrix (Pmat) used to construct the preconditioner. `Pmat` is often the same as `Amat`.
    *   Input: `ksp` - The KSP context.
    *   Input: `Amat` - The operator matrix.
    *   Input: `Pmat` - The matrix for preconditioner construction.
*   **`KSPSetUp(KSP ksp)`**: Prepares the KSP solver for solving the linear system. This includes setting up the preconditioner if one is to be used. It's often called automatically by `KSPSolve()`.
    *   Input: `ksp` - The KSP context.
*   **`KSPSolve(KSP ksp, Vec b, Vec x)`**: Solves the linear system $Ax = b$.
    *   Input: `ksp` - The KSP context.
    *   Input: `b` - The right-hand side vector.
    *   Input/Output: `x` - The initial guess vector (if `KSPSetInitialGuessNonzero()` is true); on output, it contains the solution.
*   **`KSPSolveTranspose(KSP ksp, Vec b, Vec x)`**: Solves the transposed linear system $A^T x = b$.
    *   Input: `ksp`, `b`.
    *   Input/Output: `x`.
*   **`KSPDestroy(KSP *ksp)`**: Destroys the KSP context and frees associated memory.
    *   Input/Output: `ksp` - Pointer to the KSP context.
*   **`KSPGetPC(KSP ksp, PC *pc)`**: Retrieves the preconditioner (`PC`) context associated with the KSP object. A PC is created automatically if one doesn't already exist.
    *   Input: `ksp` - The KSP context.
    *   Output: `pc` - The preconditioner context.
*   **`KSPSetTolerances(KSP ksp, PetscReal rtol, PetscReal abstol, PetscReal dtol, PetscInt maxits)`**: Sets the convergence tolerances for the KSP solver.
    *   Input: `ksp`, `rtol` (relative), `abstol` (absolute), `dtol` (divergence), `maxits` (max iterations).
*   **`KSPGetConvergedReason(KSP ksp, KSPConvergedReason *reason)`**: Gets the reason why the KSP iteration stopped (converged or diverged).
    *   Input: `ksp` - The KSP context.
    *   Output: `reason` - The `KSPConvergedReason` enum value.
*   **`KSPSetFromOptions(KSP ksp)`**: Sets KSP parameters from the PETSc options database (command-line arguments).
*   **`KSPView(KSP ksp, PetscViewer viewer)`**: Prints information about the KSP object to a viewer.
*   **`KSPGetIterationNumber(KSP ksp, PetscInt *its)`**: Gets the number of iterations performed in the last solve.
*   **`KSPSetInitialGuessNonzero(KSP ksp, PetscBool flg)`**: Informs KSP if the initial guess `x` passed to `KSPSolve()` is non-zero.
*   **`KSPSetNormType(KSP ksp, KSPNormType normtype)`**: Sets the norm type used for convergence monitoring (e.g., preconditioned, unpreconditioned, natural).
*   **`KSPSetPCSide(KSP ksp, PCSide side)`**: Sets whether the preconditioner is applied on the left, right, or symmetrically.

## Important Variables/Constants

Several enums and constants are essential for KSP configuration:

*   **`KSPType` (string enum, e.g., `KSPGMRES`, `KSPCG`)**: Defines the Krylov subspace method. Common types include:
    *   `KSPRICHARDSON`: Richardson iteration.
    *   `KSPCHEBYSHEV`: Chebyshev iteration.
    *   `KSPCG`: Conjugate Gradient method (for symmetric positive definite systems).
    *   `KSPGMRES`: Generalized Minimal Residual method.
    *   `KSPBCGS`: BiConjugate Gradient Stabilized method.
    *   `KSPTFQMR`: Transpose-Free Quasi-Minimal Residual method.
    *   `KSPPREONLY`: Applies only the preconditioner (often used with direct solvers like LU or Cholesky, accessed via the `PC` object).
    *   Many others are available.
*   **`KSPNormType` (enum `KSPNormType`)**: Specifies the norm used for convergence testing.
    *   `KSP_NORM_DEFAULT`: Default, usually `KSP_NORM_PRECONDITIONED`.
    *   `KSP_NORM_NONE`: No norm is computed explicitly by KSP (some methods track norms internally or do not require them).
    *   `KSP_NORM_PRECONDITIONED`: Norm of the preconditioned residual.
    *   `KSP_NORM_UNPRECONDITIONED`: Norm of the true residual ($b - Ax$).
    *   `KSP_NORM_NATURAL`: Natural norm, suitable for symmetric preconditioning.
*   **`PCSide` (enum `PCSide`)**: Determines how the preconditioner is applied relative to the operator.
    *   `PC_LEFT`: $P^{-1}Ax = P^{-1}b$ (default for most KSP types).
    *   `PC_RIGHT`: $A P^{-1} (Px) = b$ (default for `KSPFGMRES`).
    *   `PC_SYMMETRIC`: $P^{-1/2} A P^{-1/2} (P^{1/2}x) = P^{-1/2}b$ (e.g., for `KSPQCG`).
*   **`KSPConvergedReason` (enum `KSPConvergedReason`)**: Indicates why the KSP solver terminated. Values distinguish between various convergence criteria being met (e.g., `KSP_CONVERGED_RTOL_NORMAL`, `KSP_CONVERGED_ATOL_NORMAL`) or divergence conditions (e.g., `KSP_DIVERGED_ITS`, `KSP_DIVERGED_DTOL`, `KSP_DIVERGED_BREAKDOWN`).
*   **`KSP_CLASSID` (global `PetscClassId`)**: The internal PETSc class identifier for KSP objects, used for type checking and object management.

## Usage Examples

A conceptual C code snippet demonstrating a typical KSP usage pattern:

```c
#include <petscksp.h>

int main(int argc, char **args) {
    Mat            A;         // Operator matrix
    Vec            x, b;      // Solution and RHS vectors
    KSP            ksp;       // KSP solver context
    PC             pc;        // Preconditioner context (optional for direct manipulation)
    PetscErrorCode ierr;
    PetscInt       its;
    KSPConvergedReason reason;

    ierr = PetscInitialize(&argc, &args, (char*)0, "KSP Example"); CHKERRQ(ierr);

    // ... (Code to create and assemble Mat A, Vec b, and Vec x) ...
    // For example:
    // ierr = MatCreate(PETSC_COMM_WORLD, &A); CHKERRQ(ierr);
    // ierr = MatSetSizes(A, PETSC_DECIDE, PETSC_DECIDE, N_global, N_global); CHKERRQ(ierr);
    // ierr = MatSetFromOptions(A); CHKERRQ(ierr);
    // ierr = MatSetUp(A); CHKERRQ(ierr);
    // ... (MatSetValues and MatAssemblyBegin/End for A) ...
    // ierr = VecCreate(PETSC_COMM_WORLD, &b); CHKERRQ(ierr);
    // ierr = VecSetSizes(b, PETSC_DECIDE, N_global); CHKERRQ(ierr);
    // ierr = VecSetFromOptions(b); CHKERRQ(ierr);
    // ... (VecSetValues and VecAssemblyBegin/End for b) ...
    // ierr = VecDuplicate(b, &x); CHKERRQ(ierr); // x often starts as zero or an initial guess

    // 1. Create KSP context
    ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);

    // 2. Set operators for the linear system
    ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr); // Using A for both operator and preconditioner matrix

    // 3. (Optional) Set KSP type and PC type
    // ierr = KSPSetType(ksp, KSPGMRES); CHKERRQ(ierr);
    // ierr = KSPGetPC(ksp, &pc); CHKERRQ(ierr);
    // ierr = PCSetType(pc, PCILU); CHKERRQ(ierr);

    // 4. Set KSP options from the command line (can override KSP/PC types, tolerances, etc.)
    ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr);

    // 5. (Optional) Set specific tolerances programmatically
    // ierr = KSPSetTolerances(ksp, 1.e-6, PETSC_DEFAULT, PETSC_DEFAULT, 1000); CHKERRQ(ierr);

    // 6. Solve the linear system Ax = b
    ierr = KSPSolve(ksp, b, x); CHKERRQ(ierr);

    // 7. Check convergence and get iteration count
    ierr = KSPGetConvergedReason(ksp, &reason); CHKERRQ(ierr);
    if (reason < 0) {
        PetscPrintf(PETSC_COMM_WORLD, "KSP solver diverged for reason: %D\n", reason);
    } else {
        ierr = KSPGetIterationNumber(ksp, &its); CHKERRQ(ierr);
        PetscPrintf(PETSC_COMM_WORLD, "KSP solver converged in %D iterations.\n", its);
    }

    // 8. Destroy KSP context (this also destroys the PC by default)
    ierr = KSPDestroy(&ksp); CHKERRQ(ierr);

    // ... (Free Mat A, Vec x, Vec b) ...
    ierr = PetscFinalize();
    return 0;
}
```

## Dependencies and Interactions

KSP is a central numerical solver in PETSc and relies heavily on other components:

*   **`Mat` (Matrices)**: `KSPSetOperators()` associates one or two `Mat` objects with the KSP context. One matrix (`Amat`) defines the linear system $A$, and another (`Pmat`, often the same as `Amat`) is used for constructing the preconditioner.
*   **`Vec` (Vectors)**: KSP solvers operate on `Vec` objects, taking a right-hand side vector `b` and an initial guess/solution vector `x` as input to `KSPSolve()`.
*   **`PC` (Preconditioners)**: KSP almost always uses a `PC` object to apply preconditioning, which is crucial for the performance of iterative methods. `KSPGetPC()` allows access to the `PC` context for configuration (e.g., setting `PCType`, factorization options).
*   **`SNES` (Scalable Nonlinear Solvers)**: For solving nonlinear systems, SNES typically uses a KSP solver internally to solve the Newton linear system ($J \delta x = -F(x)$) at each nonlinear iteration. `SNESGetKSP()` provides access to this internal KSP.
*   **`TS` (Time Stepping Solvers)**: Implicit time integration methods in TS often require solving linear systems at each time step. These are usually handled by KSP, often managed internally by a SNES solver if the time-stepping problem is nonlinear.
*   **`DM` (Data Management)**: While KSP can operate on `Mat` and `Vec` objects directly, `DM`s are often used to manage the problem's underlying geometry, discretization, and parallel data distribution. `KSPSetDM()` can associate a `DM` with a KSP context, allowing the KSP (or its PC) to query the `DM` for geometric information, especially useful for multigrid preconditioners (`PCMG`).

The `src/ksp/ksp/interface/itfunc.c` file itself provides the high-level KSP interface and does not directly introduce dependencies on external (non-PETSc) libraries. Such dependencies are typically encapsulated within specific `PC` (e.g., Hypre, MUMPS via `PCLU` or `PCCHOLESKY`) or `Mat` types.
