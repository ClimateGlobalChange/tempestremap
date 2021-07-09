TempestRemap
=============

Author:  Paul Ullrich
Email:   paullrich@ucdavis.edu

Overview
--------

TempestRemap is a conservative, consistent and monotone remapping package for
arbitrary grid geometry with support for finite volumes and finite elements.
There is still quite a bit of work to be done, but any feedback is appreciated
on the software in its current form.

If you choose to use this software in your work, please cite our papers:

Paul A. Ullrich and Mark A. Taylor, 2015: Arbitrary-Order Conservative
and Consistent Remapping and a Theory of Linear Maps: Part 1.
Mon. Wea. Rev., 143, 2419–2440, doi: 10.1175/MWR-D-14-00343.1

Paul A. Ullrich, Darshi Devendran and Hans Johansen, 2016: Arbitrary-Order
Conservative and Consistent Remapping and a Theory of Linear Maps, Part 2.
Mon. Weather Rev., 144, 1529-1549, doi: 10.1175/MWR-D-15-0301.1. 

[![Install with conda](https://anaconda.org/conda-forge/tempest-remap/badges/installer/conda.svg)](https://anaconda.org/conda-forge/tempest-remap)
[![Platforms](https://anaconda.org/conda-forge/tempest-remap/badges/platforms.svg)]()

Build Instructions with Autotools
----------------------------------

TempestRemap now supports a robust autotools-based configuration/build system. Building with the standard autotools toolchain simplifies the process of linking for downstream packages like [MOAB](https://bitbucket.org/fathomteam/moab), and in [conda feedstock](https://anaconda.org/conda-forge/tempest-remap). The builds produce both linkable TempestRemap libraries and several specialized tools (or drivers) that can be installed in user-specified installation directories.

Please follow the build instructions below:

  1. `cd $TEMPESTREMAP_SRCDIR && autoreconf -fi`
  2. `mkdir -p build && cd build`
  3. Configure TempestRemap:
  ```
     ../configure --prefix=$INSTALL_DIR \ # Install dir for TempestRemap libraries
                  --with-blas=$BLAS_LIB \ # Path to BLAS libraries
                  --with-lapack=$LAPACK_LIB \ # Path to LAPACK libraries
                  --with-netcdf=$NETCDF_DIR \ # With NetCDF-C interfaces
                  --with-hdf5=$HDF5_DIR # If NetCDF was build with HDF5
  ```
  5.  Build TempestRemap: `make all`
  6.  Install TempestRemap: `make install`

Additionally, users can provide the appropriate compilers with the environmental flags (`CC`, `CXX`, `FC`, `F77`, etc.) and control the compilation/link flags (`CXXFLAGS`, `CPPFLAGS`, `LDFLAGS`, `LIBS`) as necessary.

Build Instructions with make
----------------------------

The software can be obtained from the GITHub repository via git:
```
git clone https://github.com/paullric/tempestgecore.git
```
You will likely need to edit the first couple lines of the Makefile to
customize the NetCDF paths and change any compiler flags.  Once you have
modified the Makefile, build the code:
```
make -f Makefile.gmake all
```
To clean out the object file and return the sources to pristine condition,
you can execute the following:
```
make -f Makefile.gmake clean
```

Mesh Generation
---------------

The remapping process requires multiple stages.  First you will need an Exodus
file (file extension .g) for your input mesh and your output mesh.  This can
either be done via the SQuadGen mesh utility, or via the three GenerateMesh
executables that come with TempestRemap.

For a cubed-sphere mesh:
```
./GenerateCSMesh --res <Resolution> --alt --file <Output mesh filename>.g
```
For a latitude-longitude mesh:
```
./GenerateRLLMesh --lon <longitudes> --lat <latitudes> --file <Output mesh filename>.g
```
For a geodesic mesh:
```
./GenerateICOMesh --res <Resolution> --dual --file <Output mesh filename>.g
```
Once your input and output meshes are generated, you will need to generate the
overlap mesh (that is, the mesh obtained by placing the input and output mesh
overtop one another and recalculating intersections).  This can be done as
follows:
```
./GenerateOverlapMesh --a <Input mesh>.g --b <Output mesh>.g --out <Overlap mesh>.g
```

Offline Map Generation
----------------------

Once the overlap mesh is generated, you can now generate the weight file, which
the contains information on remapping from one mesh to the other.  The type
of offline map desired is specified by `--in_type` and `--out_type`, which can be
one of the following:

fv   - Finite volume mesh, with degrees of freedom stored as volume averages
cgll - Continuous finite element method (such as spectral element)
dgll - Discontinuous finite element method (such as discontinuous Galerkin)

Offline map generation is then performed as follows:

For finite volume to finite volume remapping:
```
./GenerateOfflineMap --in_mesh <Input mesh>.g --out_mesh <Output mesh>.g \
                     --ov_mesh <Overlap mesh>.g --in_np <Remapping Order> \
                     --out_map <Output map>.nc
```
Monotone remapping in this case can be achieved with `--in_np 1`.

For finite element to finite volume remapping:
```
./GenerateOfflineMap --in_mesh <Input mesh>.g --out_mesh <Output mesh>.g \
                     --ov_mesh <Overlap mesh>.g --in_type [cgll|dgll] \
                     --out_type fv --in_np <Input order> --out_map <Output map>.nc
```
Monotone remapping in this case can be achieved with argument `--mono`.

For finite volume to finite element remapping:
```
./GenerateOfflineMap --in_mesh <Input mesh>.g --out_mesh <Output mesh>.g \
                     --ov_mesh <Overlap mesh>.g --in_type fv --out_type [cgll|dgll] \
                     --in_np <Input order> --out_np <Output order> --out_map <Output map>.nc
```
Monotone remapping in this case requires `--mono` and `--in_np 1`.

For finite element to finite element remapping:
```
./GenerateOfflineMap --in_mesh <Input mesh>.g --out_mesh <Output mesh>.g \
                     --ov_mesh <Overlap mesh>.g --in_type [cgll|dgll] \
                     --out_type [cgll|dgll] --in_np <Input order> \
                     --out_np <Output order> --out_map <Output map>.nc
```
Monotone remapping in this case requires `--in_np 1` and `--out_np 1`.

In each case, the linear weights file will then be written to `<Output map>.nc`
in SCRIP format (although it’s a bare-bones version of SCRIP format at the
moment and I’m not sure it’ll work with SCRIP utilities).  Now that the map is
generated you can apply it to your data files:

Offline Map Application
-----------------------

The offline map can be applied using the `ApplyOfflineMap` utility:
```
./ApplyOfflineMap --map <Output map>.nc --var <Comma-separated list of variables> \
                  --in_data <Input data>.nc --out_data <Output data>.nc
```
The remapped fields should then appear in `<Output data>.nc`.  Note that if your
output mesh is rectilinear, such as a latitude-longitude mesh, the data will
automatically be arranged with horizontal spatial dimensions lat and lon.

Summary
-------

Please let me know if you have any problems / bugs / comments / feature requests.


