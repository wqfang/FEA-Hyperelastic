# FEA-Hyperelastic

The latest version of elastic code is given by `NeoHookean_FinalVersion_Jul7.nb`

## Theory

I refer the theory and implementation of the neo-Hookean solid model from Prof. Bower's [Applied Mechanics of Solid](http://solidmechanics.org/Text/Chapter8_4/Chapter8_4.php) book.

## Implementation 

Some important implementation details are list below:

1. I use [`FEMAddOns`](https://github.com/WolframResearch/FEMAddOns) package to import Abaqus mesh. For installation of the package, see instructions [here](https://www.wolfram.com/language/12/nonlinear-finite-elements/contribute-fem-programmatic-utilities.html?product=language)
2. I do not use `Notation` package, which cannot be loaded on remote kernels (Brown's CCV). Therefore, it can be directly converted to a Wolfram Language Script (.wls) file and ran on CCV. The cons is that no subscriptions on symbol can be used. To run it on CCV, uncomment the first cell and comment `dir = NotebookDirectory[];  SetDirectory[dir]` in the second cell.
3. I applied **Selectively Reduced Integration** to avoid volumetric locking. Basically, volumetric part of element stiffness and residual integral is evaluated using the full integration scheme while the deviatoric part is evaluated using reduced integration points. See [here](http://solidmechanics.org/Text/Chapter8_6/Chapter8_6.php) for more details. 
4. Like the [FEA-Viscoplasticity](https://github.com/wqfang/FEA-Viscoplasticity) project, the code automatically works for both 2D and 3D geometry and hybrid types of mesh elements. 
5. `PlotMeshBC` can be used to visualize the mesh and boundary conditions.
6. I fixed a bug from the previous code. In `ApplyTraction` function, traction boundary conditions are applied on edges (faces for 3D problem). Traction components should be specified as well. Otherwise, tractions can be over-counted.

and more ...
