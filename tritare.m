function [T,TV,gamma,nodes,edges] = tritare(theta)
%function [T,TV,gamma,nodes,edges] = tritare(theta)
%
% Produces an NLEVP corresponding to planar modal vibrations
% of a Y-shaped "3-string" or "tritare"
%
% For each string, Hooke's constant (k) and linear density (rho) are 1.
% Tension in each "arm" has magnitude 1.
% Tension in the "stem" is selected to produce static equilibrium
% when the arms form the angle theta.
%
% INPUT
%     theta: angle between the "arm" strings.
%         theta should be in the interval [0,pi).
%             (default value 2*pi/3)
%         The "stem" string lies opposite the bisector of this angle
%         so that it lies on a line of symmetry.
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

if nargin < 1
    theta = 2*pi/3;
end

% Magnitude of tension in arms (strings 2 and 3)
T2 = 1;
% Magnitude of tension in stem (string 1)
T1 = 2*cos(theta/2);
% Stretched string lengths (each has relaxed length 1 and stretch factor T+1)
L1 = T1+1;
L2 = T2+1;

nodes = [  0              ,   0
         -L1              ,   0
          L2*cos(theta/2) ,  L2*sin(theta/2)
          L2*cos(theta/2) , -L2*sin(theta/2)];

edges = [2,1
         3,1
         4,1];

[T,TV,gamma] = general_web(nodes,edges);

