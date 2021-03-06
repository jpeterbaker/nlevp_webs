function [T,TV,gamma,nodes,edges] = spider2(d)
%function [T,TV,gamma,nodes,edges] = spider2(d)
%
% Orb weaver web model
% 
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem
% The network is a tracing of the web of an orb weaver spider
%
% INPUTS
%
% d is the number of spatial dimensions for the web to occupy (2 or 3)
%     default 2
%
% OUTPUTS
% 
% T is a function handle that accepts a scalar and returns a square matrix
% of dimension (2*n*ne) x (2*n*ne)
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
%     A for string 1 dimension n
%     B for string 1 dimension n
%     A for string 2 dimension 1
%     ...
%     B for string nv dimension n
%
%     a,b for each string and dimension are determined by physical parameters
%
% TV is an n x n x ne array of basis vectors
%     This is needed for understanding mode shapes
%     TV(:,i,j) is the direction of "dimension i" for string j
%
% gamma is an ne x n array of eigenvalues of the wave operator
%     gamma(i,1) is the eigenvalue for the longitudinal component of string i
%     gamma(i,j) for j=2..n is the eigenvalue for all transverse directions of string i
%          The same value is repeated n-1 times

if nargin<3 || isempty(d)
    d = 2;
end

load('orb.mat')

% If a third dimension is requested, put web in z=0 plane
if d == 3
    nodes(1,3) = 0;
end

[T,TV,gamma] = general_web(nodes,edges,d);

