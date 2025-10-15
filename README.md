# Counting spinc structures on real Bott manifolds

The applications counts the numbers of real Bott manifolds with spinc and spin structures in dimensions from **3 to 11**.

## Requirements

- OpenMP
- gcc compiler
- GNU make
- glib library
- xxhash library
- nauty library with tls support

## Compiling

Standard `configure` and `make` procedure is applicable here.

## Examples of working with RBMs

The generation of data of real Bott manifolds (RBMs), including the number of

- RBMs
- orientable RBMs
- RBMs with spin<sup>c</sup> structure
- RBMs with spin structure

in low dimensions is presented in [example](EXAMPLE.md) file.

## Applications

### Main workers of the project

- **mats**: Counts Bott matrices with spinc and spin structures, optionally prints out their d6 codes. Works in dimensions up to 10.
- **backtrack**: Counts Bott matrices with spinc and spin structures, optionally prints out their d6 codes. Works in higher dimensions also.
- **minimalf**: Filters input d6 codes for those which correspond to a minimal canonical one. This gives essentially a digraph corresponding to a diffeomorphism class of a real Bott manifold. Runs in parallel.
- **orbitg**: One-threaded version of `minimalf`, uses caching methods and hence can be very memory-consuming.

### Some helper applications

All of the following programs read from *stdin* and write to *stdout*.

- **canonicalg**: Translates list of d6 codes to the ones that represent canonical directed acyclic graphs (DAGs).
- **orientedf**: Filters input d6 codes to those which represent orientable real Bott manifolds, i.e. DAGs with all vertices of even out-degree.
- **spincf**: Filters input for matrices of spinc real Bott manifolds. *Warning:* The application assumes that the input consists of DAGs in topological order, i.e. their adjacency matrices are strictly upper triangular.
- **spinf**: Filters input for matrices of spin real Bott manifolds. *Warning:* The application assumes that the input consists of DAGs in topological order, i.e. their adjacency matrices are strictly upper triangular.
- **uniqueg**: Outputs unique elements from the input.
- **upperf**: Filter input for those d6 codes, which correspond to DAGs in topological order.
- **upperg**: Tranform every input d6 DAG code to a code of DAG in topological order. It does not check whether the input is really a DAG.

### Test applications

Test, similar to helper applications, accepts data from *stdin*.

- **test-canon**: Tests equivalence of two approaches to canonization of DAGs:

    - direct convertion of d6 codes, with output internal key128 type
    - convertion of d6 code to internal mat format, canonization and conversion to key128 type

- **test-pack**: Tests equivalence of two conversions:

    - d6 code to key128
    - d6 code to mat, followed by mat to key128

    In addition test for the mutual inversion of:

    - mat to key128
    - key128 to mat

    is taken.