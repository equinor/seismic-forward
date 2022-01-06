## GitHub

To be able to check out the repository to your file system, you must have
generated a personal access token in GitHub. This token is used in place
of your personal password when accessing the repo

  https://github.com/settings/tokens

To have Git remember the credentials for you, do this

```
git config --global credential.helper cache                  # (activate memory cache)
git config --global credential.helper 'cache --timeout=3600' # (timeout in seconds)```


## Prerequisites

  Boost     -> version 1.65.1 is included under 3rd-party
               NB! If another version of Boost is installed in the operating
               system cmake may become confused and stop working.

  Intel MKL -> Must be installed


## To make an executable using CLion or another IDE

Go to the directory containing the SeismicForward git repository

```
clion CMakeLists.txt
```


## To make an executable using cmake

Go to the directory above the directory containing the SeismicForward git repository

```
cd dir-above-SeismicForward
```
Make directory where you want the project and executable to be build

```
mkdir my-proj-dir
```

Go to this directory and run run_cmake.sh to set up compiler and library dependencies

```
cd my-proj-dir
../SeismicForward/run_cmake.sh
```

Generate the executable

```
make
```


## Run the tests

To run all tests do

```
../SeismicForward/TestScript.pl
```

To run a single test do

```
../SeismicForward/TestScript.pl case=1
```

NBNB! If more than one angle is involved, I think that results for the first
angle only is checked. This is the way NRLib currently works.
