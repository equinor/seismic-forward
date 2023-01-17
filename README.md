## Documentation

User manuals are located in the 'doc' directory.

## Prerequisites

  Boost     -> version 1.65.1 is included under 3rd-party
               NB! If another version of Boost is installed in the operating
               system cmake may become confused and stop working.

  Intel MKL -> Must be installed


## GitHub

To have Git remember the credentials for you, do this

```
git config --global credential.helper cache                  # (activate memory cache)
git config --global credential.helper 'cache --timeout=3600' # (timeout in seconds)```
```

## Fork and checkout a branch

If you intend to modify the repository, you should checkout a fork of this
project rather than checking out the repository directly. When you have made
a fork, you can set this repository as your origin using:

```
git checkout --track origin/release-4.3
```
(At the time being, it is adviced that you use branch release-4.3)


## To make an executable using CLion or another IDE

Go to the directory containing the seismic-forward git repository

```
clion CMakeLists.txt
```


## To make an executable using cmake

Go to the directory above the directory containing the seismic-forward git repository

```
cd dir-above-seismic-forward
```
Make directory where you want the project and executable to be build

```
mkdir my-proj-dir
```

Go to this directory and run run_cmake.sh to set up compiler and library dependencies

```
cd my-proj-dir
../seismic-forward/run_cmake.sh
```

Generate the executable

```
make
```


## Run seismic-forward

The executable can be run directly once an XML model file/job file has been created. The only input data needed is an Eclipse-grid with Vp, Vs and Rho.

```
cd directory-with-Eclipse-grid
path-to-executable/seismic_forward my_job_file.xml
```


## Run the tests

To run all tests do

```
../seismic-forward/TestScript.pl
```

To run a single test do

```
../seismic-forward/TestScript.pl case=1
```

NBNB! If more than one angle is involved, I think that results for the first angle only is checked. This is the way NRLib currently works.
