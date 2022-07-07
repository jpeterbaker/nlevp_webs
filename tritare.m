function [T,TV,gamma,nodes,edges,eigs] = tritare()
%function [T,TV,gamma,nodes,edges,eigs] = tritare()
%
% Produces an NLEVP corresponding to planar modal vibrations of a 3-string or "tritare"
%
% For this problem
%    Nodes        nv = 4
%    Strings      ne = 3
%    Dimensions   d  = 2
%
% OUTPUTS
% 
% T is a function handle that accepts a scalar and returns a square matrix
% of dimension (2*d*ne) x (2*d*ne)
%     If T(omega) is singular, then omega/(2*pi) is a natural frequency of the web
%     The corresponding singular vector gives coefficients of mode shape
%     The mode shape is described as a function X(x) for each string
%     Y(x) is X(x) transformed to a coordinate system parallel to the string
%     Each Y(x) has the form A*sin(a*x)+B*cos(b*x) in each spatial dimension
%
%     The coefficient order is
%     A for string 1 dimension 1
%     B for string 1 dimension 1
%     A for string 1 dimension 2
%     ...
%     A for string 1 dimension d
%     B for string 1 dimension d
%     A for string 2 dimension 1
%     ...
%     B for string nv dimension d
%
%     a,b for each string and dimension are determined by physical parameters
%
% TV is an d x d x ne array of basis vectors
%     This is needed for understanding mode shapes
%     TV(:,i,j) is the direction of "dimension i" for string j
%
% gamma is an ne x d array of eigenvalues of the wave operator
%     gamma(i,1) is the eigenvalue for the longitudinal component of string i
%     gamma(i,j) for j=2..d is the eigenvalue for all transverse directions of string i
%          The same value is repeated d-1 times
%
% nodes is an nv x d matrix with the coordinates of nodes
%     nodes(i,:) is the location of node i
%     nv is the number of nodes
%     d is the number of spatial dimensions (usually 2 or 3)
% 
% edges is an ne x 2 matrix with the indices of connected nodes
%     edges(i,:) contains the indices of the two nodes connected by string i
%     ne is the number of edges
%
%
% eigs is a 24-vector of the smallest positive eigenvalues
%
s2 = 1/sqrt(2);
nodes = [   0,  0
           -1,  0
           s2, s2
           s2,-s2];

edges = [2,1
         3,1
         4,1];

[T,TV,gamma] = general_web(nodes,edges);

% Natural frequencies are solutions of
%
% (1)   sin(x) = 0
%
% (2)   sin(x/sqrt(2)) = 0
%
% (3)   tan(x) + sqrt(2)*tan(x/sqrt(2)) = 0
%    or equivalently
%       sin(x)*cos(x/sqrt(2)) + sqrt(2)*sin(x/sqrt(2))*cos(x) = 0
%
% (4)   2*sqrt(2)*tan(x) + tan(x/sqrt(2)) = 0
%    or equivalently
%       2*sqrt(2)*sin(x)*cos(x/sqrt(2)) + cos(x)*sin(x/sqrt(2)) = 0
% 
% These solutions are not difficult to find with numerical root-finding methods

eigs = [1.78993037359553312790910971176
        1.99658732056910586526458512505
        3.14159265358979323846264338328
        3.43723396871777041008382193164
        3.77073144615781997855778892843
        4.44288293815836624701588099006
        5.41991645432706565399912106673
        5.80985767635318626426947322806
        6.28318530717958647692528676656
        7.10959048557161358103777053585
        7.44389306692523896182785166003
        8.88576587631673249403176198012
        9.15688541405701741654222854053
        9.31460464998877656090308499856
        9.42477796076937971538793014984
        11.0327967546384344713384388967
        11.0699562440765798653039052732
        12.5663706143591729538505735331
        12.7256215934544782848147064046
        12.9429962493987284344002807234
        13.3286488144750987410476429702
        14.6305702085998011895379732584
        14.9980859062085563610159486292
        15.7079632679489661923132169164];
