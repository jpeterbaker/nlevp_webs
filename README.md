# Nonlinear Eigenvalue Problems for vibrating webs

> Can a spider tune its web to sing in the wind?

The natural frequencies and mode shapes of multibody structures can be found
by collecting the connection constraints into a nonlinear eigenvalue problem (NLEVP).
This project provides tools for constructing such NLEVPs
for networks of elastic strings such as spiderwebs.
The NLEVPs in this project are similar to dynamic stiffness matrices,
but they do not need to be symmetric
and their derivation does not require algebra to eliminate displacement constraints.

## Network format

Networks of strings are described by two variables: `nodes` and `edges`.
* `nodes` is double-precision matrix of dimension ne-by-d where
  * ne is the number of strings (edges)
  * d is the number of spatial dimensions (usually 2 or 3)
* `edges` is an integer matrix of node connections of dimension ne-by-2
  * Row `i` has two indices indicating two nodes that are connected by a string
  * The indices in `edges` refer to the rows of `nodes`

## Files

The code files of the project are listed below with brief descriptions.

### Main function

* `general_web.m`

This is the most powerful function in the project.
It constructs the NLEVP of arbitrary string networks.

* `incidence_vectors.m`

This is a utility function for preprocessing the network description variables `nodes` and `edges`.

### Tritare

* `tritare.m`

Produces the NLEVP of a 3-string or "tritare."

* `mode_curves.m`

Produces curved string shapes during modal vibration.
If the result does not appear to be continuous,
then the eigenvalue and eigenvector inputs are probably not solutions
to the NLEVP.

* `demo_tritare.m`

Demonstrates how the outputs of `tritare.m` and `mode_curves.m`
may be used to draw a modally-vibrating 3-string.

### Regular webs

* `regweb.m`

Produces the NLEVP of a web of regularly-spaced spokes and rings.
The number of spokes and rings are customizable,
but 5 spokes and 4 rings is the default.

* `regweb_graph.m`

This is a utility function that produces the 

### Realistic webs

* `spider1.m`

Provides the NLEVP and network description of a Deinopis web.

* `spider2.m`

Provides the NLEVP and network description of an orb-weaver web.

* `deinopis.mat`

Data file with a tracing of a Deinopis web.

* `orb.mat`

Data file with a tracing of an orb-weaver web.

