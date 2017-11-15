# Metos3D

**Marine Ecosystem Toolkit for Optimization and Simulation in 3-D**

## Installer

```
$>
curl -O https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh

conda install -c jpicau metos3d
```

## Base system

You've cloned the base Metos3D system. This includes ...

```
metos3d         # work horse, Python script
src/            # transport driver, C code
model/          # bgc models, Fortran code
```

## Environment

Metos3D needs a sane compile environment, this includes ...

*If you are working on a known machine and want to use one of the standard models,
have a look at our [environment section](), we have prepared a set of standard versions,*

however, if you want the full control, you'll need,
*sequential* execution:
C/C++, Fortran compiler,
*parallel* execution,
MPI versions of C/C++,

with those compilers,
a compiled PETSc library,
version 3.7., which means,
two environment variables set,
named `PETSC_DIR` and `PETSC_ARCH`,

## Data

You're missing the data, yet. You can clone it.
However, you'll need Git LFS (large file support).

```
$> git clone http://
```

or
you can download an compressed archive, using `curl` or `wget`

```
$> wget https://
$> curl -O https://
```
unzip,


## Compile and run

Now, you can ...
compile an executable and run a simulation

for instance:

```
$> metos3d simpack NPZD-DOP
$> ./metos3d-simpack-NPZD-DOP.exe model/NPZD-DOP/option/test.NPZD-DOP.option.yaml
```

see [Metos3D cheat sheet]() for a detailed reference of the subcommands of the `metos3d` script,




# Introducing Metos3D version 1.0.0
following [sematic versioning]() (**major.minor.patch**),
change in major version means loss of compatibility

we hybernated the development of Metos3D in separated repositories,
the last avaiable versions were: `metos3d v0.6.0`, `simpack v0.5.0`, `model v0.4.0` and `data v0.2.1`
all git repositories were archieved as one file and can be downloaded [here](),
each repo includes its earlier versions,
in particular access to versions used in [Piwonski and Slawig, 2016]()

# New concept
Metos3D 1.0.0 introduces a new concept of usage,
one repository,

we distinguish two installation types:

- container, easy to install and use, hard-coded,
    one executable for each model, inflexible,
    including everything, huge 
- script/system, requires [mpi versions of c and fortran compilers](),
    as well as a compiled version of [petsc](),
    fully flexible, lean

# Documentation
- cheat sheet, BGC cheat sheet, templates
- tutorial
- manual
- code reference

# References
- [Piwonski and Slawig, 2016]() Metos3D: ...




<!---->
<!--DUMP-->
<!---->

<!--<table>-->
<!--<tr>-->
<!--<td>metos3d</td><td>v0.6.*</td>-->
<!--</tr>-->
<!--<tr>-->
<!--<td>simpack</td><td>v0.6.*</td>-->
<!--</tr>-->
<!--<tr>-->
<!--<td>model  </td><td>v0.4.*</td>-->
<!--</tr>-->
<!--<tr>-->
<!--<td>data   </td><td>v0.2.*</td>-->
<!--</tr>-->
<!--</table>-->
~~~~
