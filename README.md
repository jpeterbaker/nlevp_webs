# Nonlinear Eigenvalue Problems for vibrating webs

<div style="float:right;margin:10px 10px 10px 30px;font-size:12px;position:relative">
The spider as an artist
<br>
&nbsp;&nbsp;&nbsp;&nbsp;Has never been employed
<br>
Though his surpassing merit
<br>
&nbsp;&nbsp;&nbsp;&nbsp;Is freely certified
<br>
<i style="right:0px;position:absolute">Emily Dickinson</i>
<br>
</div>

The natural frequencies and mode shapes of multibody structures can be found
by collecting the connection constraints into a nonlinear eigenvalue problem (NLEVP).
This project provides tools for constructing such NLEVPs
for networks of elastic strings, spiderwebs being typical examples.
The NLEVPs in this project are similar to dynamic stiffness matrices,
but they do not need to be symmetric
and their derivation does not require algebra to eliminate displacement constraints.

An NLEVP maps complex scalars to square matrices.
When the matrix is singular, the input is called an eigenvalue,
and all corresponding null vectors of the matrix are its eigenvectors.
The NLEVPs produced by this project have eigenvalues at natural frequencies of the related multibody system.
With some work, the corresponding mode shapes can be derived from the eigenvectors.

The primary purpose of this project is to provide interesting problem examples,
*not* to provide solution methods.

## Getting started

1. Download the project and update MatLab's working directory (or search path)
2. Run `demo_tritare`
   * Figure 1 should display a wiggly Y-shape that is a modal vibration of a "tritare" string
3. Examine the code of `demo_tritare.m` to see how the outputs of NLEVP-producing functions like `tritare` are used
4. Modify `demo_tritare.m` to choose a different mode (choose a value of `j` from 1 to 24)
5. Rerun `demo_tritare` to see a different mode shape

`tritare.m` is the only file in the project that provides pre-calculated eigenvaleus.
In order to draw the mode shapes of the other systems, you will need to calculate the eigenvalues yourself.

## Requirements

This project runs in base MatLab without any additional packages needed.
It was primarily tested with version 2020a.
