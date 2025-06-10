# DM (Data Management)

## Overview

The `DM` (Distribution Manager or Data Management) object is a fundamental component within PETSc designed to manage the complexities of data distribution and communication required for parallel numerical computations. It acts as an abstraction layer, mediating between the geometric and topological description of a problem's domain (like meshes or grids) and the algebraic objects (vectors and matrices) used by solvers. `DM`s are central to organizing data for various parallel data structures, enabling PETSc's solvers (KSP, SNES, TS) to work efficiently with structured and unstructured grids/meshes.

## Key Components (Functions)

The `DM` interface provides a rich set of functions for creating, configuring, and using data management objects. Here are some of the major public functions found primarily in `src/dm/interface/dm.c` and its associated header `petscdm.h`:

*   **`DMCreate(MPI_Comm comm, DM *dm)`**: Creates an empty `DM` object associated with a given MPI communicator. The specific type of `DM` (e.g., `DMDA` for structured grids, `DMPLEX` for unstructured meshes) must be set subsequently using `DMSetType()`.
    *   Input: `comm` - The MPI communicator for the `DM` object.
    *   Output: `dm` - The newly created `DM` object.
*   **`DMSetType(DM dm, DMType method)`**: Sets the specific implementation type for the `DM` object.
    *   Input: `dm` - The `DM` object.
    *   Input: `method` - A string representing the `DMType` (e.g., `DMDA`, `DMPLEX`).
*   **`DMSetUp(DM dm)`**: Performs all the setup operations for the `DM` data structures. This is typically called after various parameters and options (like dimension, distribution, stencil type for `DMDA`, or topology for `DMPLEX`) have been set.
    *   Input: `dm` - The `DM` object.
*   **`DMDestroy(DM *dm)`**: Destroys the `DM` object and frees all associated data.
    *   Input/Output: `dm` - Pointer to the `DM` object to be destroyed.
*   **`DMView(DM dm, PetscViewer viewer)`**: Visualizes or prints information about the `DM` object using a PETSc viewer.
    *   Input: `dm` - The `DM` object.
    *   Input: `viewer` - The `PetscViewer` (e.g., `PETSC_VIEWER_STDOUT_WORLD`).
*   **`DMCreateGlobalVector(DM dm, Vec *vec)`**: Creates a global vector whose parallel layout is compatible with the `DM`. Global vectors represent the distributed solution or right-hand side data without ghost points.
    *   Input: `dm` - The `DM` object.
    *   Output: `vec` - The created global `Vec`.
*   **`DMCreateLocalVector(DM dm, Vec *vec)`**: Creates a local vector. Local vectors typically include space for ghost points necessary for local computations (e.g., stencil operations).
    *   Input: `dm` - The `DM` object.
    *   Output: `vec` - The created local `Vec`.
*   **`DMGetLocalToGlobalMapping(DM dm, ISLocalToGlobalMapping *ltog)`**: Retrieves the mapping from local (including ghost points) to global indexing for the `DM`.
    *   Input: `dm` - The `DM` object.
    *   Output: `ltog` - The local-to-global mapping object.
*   **`DMRefine(DM dm, MPI_Comm comm, DM *dmf)`**: Refines the `DM` (mesh/grid) to create a finer representation suitable for multigrid methods or adaptive refinement.
    *   Input: `dm` - The `DM` to refine.
    *   Input: `comm` - Communicator for the new `DM` (or `MPI_COMM_NULL`).
    *   Output: `dmf` - The refined `DM`.
*   **`DMCoarsen(DM dm, MPI_Comm comm, DM *dmc)`**: Coarsens the `DM` to create a coarser representation.
    *   Input: `dm` - The `DM` to coarsen.
    *   Input: `comm` - Communicator for the new `DM` (or `MPI_COMM_NULL`).
    *   Output: `dmc` - The coarsened `DM`.
*   **`DMSetFromOptions(DM dm)`**: Sets parameters in the `DM` from the PETSc options database (command-line arguments).
*   **`DMCreateMatrix(DM dm, Mat *mat)`**: Creates a matrix whose parallel layout corresponds to the `DM`, suitable for operators like Jacobians.
    *   Input: `dm` - The `DM` object.
    *   Output: `mat` - The created `Mat`.
*   **`DMGetLocalSection(DM dm, PetscSection *section)`**: Gets the `PetscSection` that defines the layout of local data (degrees of freedom) over the mesh points.
*   **`DMSetLocalSection(DM dm, PetscSection section)`**: Sets the `PetscSection` for local data layout.
*   **`DMGetGlobalSection(DM dm, PetscSection *section)`**: Gets the `PetscSection` for global data layout.
*   **`DMGlobalToLocalBegin(DM dm, Vec g, InsertMode mode, Vec l)` / `DMGlobalToLocalEnd(DM dm, Vec g, InsertMode mode, Vec l)`**: Initiates and completes the scatter operation from a global vector `g` to a local vector `l` (including ghost points).
*   **`DMLocalToGlobalBegin(DM dm, Vec l, InsertMode mode, Vec g)` / `DMLocalToGlobalEnd(DM dm, Vec l, InsertMode mode, Vec g)`**: Initiates and completes the gather/accumulation operation from a local vector `l` to a global vector `g`.
*   **`DMGetCoordinatesLocal(DM dm, Vec *c)`**: Gets a local vector containing the coordinates of the mesh/grid points.
*   **`DMGetCoordinateDim(DM dm, PetscInt *dim)`**: Gets the dimension of the coordinates.

## Important Variables/Constants

The `DM` interface utilizes several important enumerated types and constants:

*   **`DMType` (string enum defined with `DMRegister`)**: Specifies the underlying implementation of the `DM`. Key types include:
    *   `DMDA`: For structured, Cartesian-like grids.
    *   `DMPLEX`: For unstructured meshes, supporting complex topologies.
    *   `DMCOMPOSITE`: For combining multiple `DM` objects.
    *   `DMSWARM`: For particle methods.
    *   `DMNETWORK`: For data on networks/graphs.
*   **`DMBoundaryType` (enum, string array `DMBoundaryTypes[]` in `dm.c`)**: Describes boundary conditions for `DMDA` and other types, e.g., `DM_BOUNDARY_NONE`, `DM_BOUNDARY_GHOSTED`, `DM_BOUNDARY_PERIODIC`, `DM_BOUNDARY_MIRROR`.
*   **`DMAdaptationType` (enum `DMAdaptationType` in `petscdm.h`)**: Indicates the type of mesh adaptation, e.g., `DM_ADAPTATION_NONE`, `DM_ADAPTATION_REFINE`.
*   **`DMPointLocationType` (enum `DMPointLocationType` in `petscdm.h`)**: Used for strategies when locating points within a `DM`.
*   **`DMPolytopeType` (enum, string array `DMPolytopeTypes[]` in `dm.c`)**: Defines various cell/element shapes, e.g., `DM_POLYTOPE_TRIANGLE`, `DM_POLYTOPE_QUADRILATERAL`, `DM_POLYTOPE_TETRAHEDRON`, `DM_POLYTOPE_HEXAHEDRON`.
*   **`DM_CLASSID` (global `PetscClassId`)**: A class identifier used internally by PETSc to identify `DM` objects. This is fundamental for PETSc's object-oriented design and type checking.
*   **`DMLABEL_CLASSID` (global `PetscClassId`)**: Class identifier for `DMLabel` objects, which are used with `DM`s to mark regions, boundaries, or subsets of mesh entities.

## Usage Examples

A conceptual C code snippet illustrating a typical `DM` lifecycle:

```c
#include <petscdm.h>
#include <petscdmplex.h> // Example for DMPLEX
#include <petscksp.h>    // For KSP (example solver)

int main(int argc, char **argv) {
    DM             dm;
    Vec            globalVec, localVec;
    Mat            A;
    // KSP            ksp; // Example solver
    PetscErrorCode ierr;
    PetscInt       dim = 2; // Example dimension

    ierr = PetscInitialize(&argc, &argv, (char*)0, "DM Usage Example"); CHKERRQ(ierr);

    // 1. Create a DM object
    ierr = DMCreate(PETSC_COMM_WORLD, &dm); CHKERRQ(ierr);

    // 2. Set the DM type and configure it
    // Example using DMPLEX for a simple box mesh:
    ierr = DMSetType(dm, DMPLEX); CHKERRQ(ierr);
    // Specific setup for DMPLEX (or DMDA, etc.)
    // ierr = DMPlexCreateBoxMesh(PETSC_COMM_WORLD, dim, PETSC_FALSE, NULL, NULL, NULL, NULL, &dm); CHKERRQ(ierr);
    // or more detailed setup:
    ierr = DMSetDimension(dm, dim); CHKERRQ(ierr);
    // ... (other DMPlex specific setup like DMPlexSetChart, DMPlexSymmetrize, etc.)


    // 3. Set options from the command line (can override type and other settings)
    ierr = DMSetFromOptions(dm); CHKERRQ(ierr);

    // 4. Finalize DM setup after all options are set
    ierr = DMSetUp(dm); CHKERRQ(ierr);

    // (Optional) Distribute the DM if it's a DMPLEX and parallel
    // DM dmDist = NULL;
    // DMPlexDistribute(dm, 0, NULL, &dmDist);
    // if (dmDist) {
    //   DMDestroy(&dm);
    //   dm = dmDist;
    // }

    // 5. Create global and local vectors associated with the DM
    ierr = DMCreateGlobalVector(dm, &globalVec); CHKERRQ(ierr);
    ierr = DMCreateLocalVector(dm, &localVec); CHKERRQ(ierr);

    // 6. Create a matrix associated with the DM (e.g., for a PDE operator)
    ierr = DMCreateMatrix(dm, &A); CHKERRQ(ierr);

    // 7. Use the DM with solvers (e.g., KSP, SNES, TS)
    //    The DM is often passed to the solver context:
    // ierr = KSPCreate(PETSC_COMM_WORLD, &ksp); CHKERRQ(ierr);
    // ierr = KSPSetDM(ksp, dm); CHKERRQ(ierr);
    // ierr = KSPSetOperators(ksp, A, A); CHKERRQ(ierr);
    // ... (setup KSP further and solve) ...

    // 8. Perform data exchange operations
    // ierr = DMGlobalToLocalBegin(dm, globalVec, INSERT_VALUES, localVec); CHKERRQ(ierr);
    // ierr = DMGlobalToLocalEnd(dm, globalVec, INSERT_VALUES, localVec); CHKERRQ(ierr);
    // ... (work with localVec, which now has updated ghost point data) ...


    // 9. Clean up
    ierr = MatDestroy(&A); CHKERRQ(ierr);
    ierr = VecDestroy(&localVec); CHKERRQ(ierr);
    ierr = VecDestroy(&globalVec); CHKERRQ(ierr);
    // ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
    ierr = DMDestroy(&dm); CHKERRQ(ierr);

    ierr = PetscFinalize();
    return ierr;
}
```

## Dependencies and Interactions

The `DM` subsystem is a core part of PETSc and has strong dependencies and interactions with several other components:

*   **`Vec` (Vectors)**: `DM`s are responsible for creating `Vec` objects (`DMCreateGlobalVector`, `DMCreateLocalVector`) that have parallel layouts consistent with the `DM`'s decomposition of the domain. `VecGetDM()` and `VecSetDM()` provide a link between vectors and their managing `DM`.
*   **`Mat` (Matrices)**: `DM`s define the parallel structure for matrices (`DMCreateMatrix`), which are typically used to represent linear operators (e.g., Jacobians in SNES, discretization matrices in KSP). `MatGetDM()` and `MatSetDM()` link matrices to their `DM`.
*   **`PetscSF` (StarForest)**: `PetscSF` is used by `DM` for managing communication patterns, especially for updating ghost point data in vectors (`DMGetPointSF`, `DMGetSectionSF`) and for mesh redistribution.
*   **`PetscSection`**: `PetscSection` objects are crucial for defining the layout of data (degrees of freedom) over the mesh points or cells. `DM`s use `PetscSection`s (`DMGetLocalSection`, `DMGetGlobalSection`) to manage how field data is distributed and numbered.
*   **`DMLabel`**: Labels are used to mark subsets of the mesh or grid (e.g., boundaries, different material regions). `DM`s manage a list of labels and use them for various operations, including defining regions for boundary conditions or specific physics.
*   **`PetscDS` (Discrete System)**: `PetscDS` objects store information about the (system of) partial differential equations being solved on the domain represented by the `DM`, including information about fields, discretizations, and boundary conditions. `DMGetDS()` and `DMCreateDS()` link `DM`s to `PetscDS` objects.
*   **Solvers (KSP, SNES, TS)**:
    *   **KSP (Linear Solvers)**: `DM`s provide the structure for matrices and vectors used by KSP. Critically, `DM`s are essential for geometric multigrid preconditioners (`PCMG`), where the hierarchy of grids is managed by a hierarchy of `DM` objects.
    *   **SNES (Nonlinear Solvers)**: `DM`s provide the structure for Jacobian matrices and residual vectors. SNES uses `DM`s for problem definition and, like KSP, for geometric multigrid methods (e.g., FAS).
    *   **TS (Timesteppers)**: For time-dependent problems, `DM`s manage the data at each time step and can be used by TS to solve the underlying (often nonlinear) systems.

External library dependencies are generally indirect for the core `dm.c` file. Specific `DM` implementations (like `DMMOAB` for MOAB or `DMPLEX` interfacing with Gmsh/Triangle via file I/O) might introduce such dependencies. LibCEED is an optional dependency for performant finite element calculations, and `DM` can manage the LibCEED context (`DMGetCeed`).
