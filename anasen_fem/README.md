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
* v0.0.2, planned TODO
	- Garfield to take Elmer results and perform charge-transport
