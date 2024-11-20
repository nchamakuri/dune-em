Preparing the Sources
=========================

 You will need to have the following programs installed on your system: 

  `cmake >= 3.13`  
  `gcc >= 9`  

The module also would require additional softwares:  
+ MPI (OpenMPI or MPICH)  
+ SuiteSparse
+ Boost libraries for C++
+ SuperLU
+ BLAS, LAPACK packages
+ HDF5 library
  
First, download the DUNE core modules to a single directory on your computer.
For example, create a directory and clone the core modules and extension modules required for `dune-heart` into it:  
_Note: Checkout release branch 2.9 when the modules are downloaded directly from git._  
```bash
mkdir Dune
cd Dune
git clone https://gitlab.dune-project.org/core/dune-common.git  
git clone https://gitlab.dune-project.org/core/dune-geometry.git  
git clone https://gitlab.dune-project.org/core/dune-localfunctions.git  
git clone https://gitlab.dune-project.org/core/dune-istl.git  
git clone https://gitlab.dune-project.org/core/dune-grid.git  
git clone https://gitlab.dune-project.org/staging/dune-uggrid.git  
git clone https://gitlab.dune-project.org/staging/dune-functions.git  
git clone https://gitlab.dune-project.org/staging/dune-typetree.git  
git clone https://gitlab.dune-project.org/extensions/dune-alugrid.git  
git clone https://gitlab.dune-project.org/pdelab/dune-pdelab.git
```

Getting started
---------------
Download the `dune-heart` module into the same directory (`Dune`) where the other modules are located.
If all the prerequisites are met, run the following command:

```bash    
     ./dune-common/bin/dunecontrol all
```

which will find all installed dune modules as well as all dune modules
(not installed) which sources reside in a subdirectory of the current
directory. Note that if Dune is not installed properly you will either
have to add the directory where the `dunecontrol` script resides (probably
./dune-common/bin) to your path or specify the relative path of the script.

Most probably, you'll have to provide additional information to dunecontrol
(e. g. compilers, configure options) and/or make options.

The most convenient way is to use options files in this case. The files
define four variables:

CMAKE_FLAGS      flags passed to cmake (during configure)  
CONFIGURE_FLAGS  

An example options file might look like this:

#use this options to configure and make if no other options are given  
CMAKE_FLAGS=" \
-DCMAKE_CXX_COMPILER=g++-9 \
-DCMAKE_CXX_FLAGS='-Wall -pedantic' \
-DCMAKE_INSTALL_PREFIX=/install/path" #Force g++-9 and set compiler flags

If you save this information into `example.opts` you can pass the opts file to
`dunecontrol` via the --opts option, e. g.
```bash
  ./dune-common/bin/dunecontrol --opts=example.opts all
```
Compling the program:
--------------------
Navigate to the build directory inside the `dune-heart` module and then to the electromechanics directory to run the cardiac mechanics simulation as follows:  
```bash
cd dune-heart/build-cmake/dune/heart/electromechanics/
```

Go to the `em_partitioned` directory for the partitioned scheme or the `em_coupled` directory for the monolithic/coupled scheme.
  ```bash
cd em_partitioned
(or)
cd em_coupled

make

./electromechanical paramL2space/electromechanical2D_80.ini
(or)
./electromechanics paramL2space/electromechanical2D_80.ini
```
To run parallelly with MPI, 
```bash
mpirun -np <number_of_processors> ./electromechanical paramL2space/electromechanical2D_80.ini
```
The `<number_of_processors>` has to be equal to the product of `xprocs` and `yprocs` in the parameter file `electromechanical2D_80.ini`.

More info
---------

See

     dunecontrol --help

for further options.  
The full build system is described in the dune-common/doc/buildsystem (Git version) or under share/doc/dune-common/buildsystem if you installed DUNE!


Contact
-------
For questions, issues, or further support, please reach out to:  
- Dr. Nagaiah Chamakuri (nagaiah.chamakuri@iisertvm.ac.in)
- Gopika P B (gopikapb23@iisertvm.ac.in)

  
