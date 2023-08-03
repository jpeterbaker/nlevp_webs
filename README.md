# Nonlinear Eigenvalue Problems for vibrating webs

<!--
Security restrictions on GitHub required me to use some some creativity
to get this poem formatted as I wanted in the preview.
-->
<div align="right">
<table>
<tr>
<td>
<pre>The spider as an artist
    Has never been employed
Though his surpassing merit
    Is freely certified
            <i>Emily Dickinson</i>
</pre>
</td>
</tr>
<table>
</div>

<div align="left">

The natural frequencies and mode shapes of multibody structures can be found
by collecting the connection constraints into a nonlinear eigenvalue problem (NLEVP).
This project provides tools for constructing such NLEVPs
for networks of elastic strings, spiderwebs being typical examples.
The NLEVPs in this project are similar to dynamic stiffness matrices,
but they do not need to be symmetric (damping is allowed)
and their derivation does not require algebra to eliminate displacement constraints.

An NLEVP maps complex scalars to square matrices.
When the matrix is singular, the input is called an eigenvalue,
and all corresponding null vectors of the matrix are its eigenvectors.
The NLEVPs produced by this project have eigenvalues at natural frequencies of the related multibody system.
The mode shapes of the structure can then be derived from the eigenvectors.

The primary purpose of this project is to provide interesting problem examples,
*not* to provide solution methods or provide realistic
natural frequencies of actual spiderwebs.

## Getting started

### Setup

1. Download the project 
2. Update MatLab's working directory (or search path) to include the project directory
3. Run `demo_tritare` to test your setup
   * A figure should display a wiggly Y-shape that is a modal vibration of a "tritare" string

### First experiment

1. Modify `demo_tritare.m` to choose a different mode (choose a value of `j` from 1 to 24)
2. Rerun `demo_tritare` to see the corresponding mode shape

You may want to examine the code of `demo_tritare.m`
and the documentation of `mode_curves.m` to see how the outputs
of NLEVP-producing functions like `tritare` can be  used.

## More reading

Notes further explaining this codebase are in preparation.
A link will be provided when ready.

## Requirements

This project runs in base MatLab without any additional packages needed.
It was developed with versions 2020a and 2021b
but is expected to work in much older versions.

</div>
