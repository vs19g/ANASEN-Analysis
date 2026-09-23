### README for ANASEN fem simulations:

* There are a few iterations of these simulations that already exist. Be sure to also locate and refer to them if necessary.
* Install gmsh and its python api by running (Ubuntu 22.04 LTS)

```
		sudo apt install gmsh python3-gmsh 
```

* Gmsh gives us the tools to create a meshgrid that samples the 2d space appropriately to plot the field/equipotential lines.
* The output file typically has the .msh extension. This is read as input to Elmer, which is the FEM differential-equation solver.

* Install Elmer via the following steps:

```
		sudo add-apt-repository ppa:elmer-csc-ubuntu/elmer-csc-ppa
		sudo apt install elmerfem-csc-eg
```

* Install ParaView for visualizations by downloading from the Linux .tar.gz link at https://www.paraview.org/download/
	- The current version is tested to work on Paraview 6.1.0. The default version in Ubuntu 22.04 repositories has some trouble with scripting
* v0.0.1, March 10 2026
	- 2d simulations of fields only. gmsh for meshing, elmer for fem, paraview to plot
	- Before running, open `paraview_plotter.py` to make the bash shebang (#!) point to the location of `pvpython` or `pvbatch`
	- `python3 run.py` should run everything in order, and is hopefully all the files are self-documenting
* v0.0.2, September 2026
	- Adds the Ramo weighting field and the `E . E_w` dot product for one cathode wire.
	- `wires_gmsh2d_bc.py` now cuts the wires into the gas disk with `occ.fragment` instead of `mesh.embed`, so each wire boundary is a real edge of the domain. One cathode (optional 2nd argument, default 1) gets its own physical group, **tag 40**; the other 23 stay in tag 30.
	- Two solves, same mesh, identical apart from the potentials:

	| file | tag 40 | output |
	|---|---|---|
	| `wires2d.sif` | 0 V, like every other cathode | `wires2d/elfield_anasen_t0001.vtu` |
	| `wires2d_weight.sif` | 1 V, all else 0 V | `wires2d/elfield_weight_t0001.vtu` |

	- `dotproduct.py` reads both and writes `wires2d/dotproduct.vtu` with `DotProduct = E . E_w` (V/m^2). It aborts if the node sets differ or if either field is all zeros — both are silent failures otherwise.
	- `paraview_dotproduct.py` renders it. `DOT_MIN`/`DOT_MAX` set the colour range and decide whether the picture shows anything.
	- `paraview_plotter.py` gains two views: the equipotentials over one quadrant, and streamlines over the same quadrant. Streamlines are drift paths up to diffusion, so they show which wire collects charge from where; `SEED_POTENTIAL` and `SEED_STRIDE` control them.
	- Six PNGs per z-locus, archived as before (`<Stem>_z_<count>_<z>[_quarter].png`). `Field_ouput` keeps its typo so the new files sort with the dozen already in `png/`.
	- `STAGE` in `run.py`: `1` = mesh + physical solve, `2` = weighting solve + dot product reusing stage 1's mesh, `0` = both. Stage 2 never re-meshes, so it applies to whichever z stage 1 ran last.
	- **Do not pass `-autoclean` to ElmerGrid.** It renumbers the physical groups (13/10/20/30/40 become 1/4/5/6/7), so every `Target Bodies`/`Target Boundaries` in the sifs matches nothing, Elmer solves nothing, and the potential comes out identically zero with no error.
	- `mesh.recombine()` and `mesh.refine()` are both off: recombination ran over half an hour on the barrel surface without finishing, and refine took the mesh to ~16M nodes.
	- Note on the solver stack: both sifs are the same file bar the potentials, so whatever `FluxSolver` does to `Electric Field` it does equally to both and the dot product stays consistent.
* v0.0.3, planned TODO
	- Garfield to take Elmer results and perform charge-transport
	- Sweep the weighting solve over all 24 cathodes rather than one at a time
