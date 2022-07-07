function [T,TV,gamma] = regweb_6_4(d)
%function [T,TV,gamma] = regweb_6_4(d)
%
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem
% The network is a "spider web" with 6 spokes and 4 rings
%
% For this problem
%    Nodes        nv = 31
%    Strings      ne = 54
%
% INPUTS
%
% d is the number of spatial dimensions for the web to occupy (2 or 3)
%     default 2
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
if nargin<1 || isempty(d)
    d=2;
end
[T,TV,gamma] = regweb(6,4,d);


