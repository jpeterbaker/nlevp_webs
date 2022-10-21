
% This script demonstrates use of a basic integration solver

% The NLEVP function
T = tritare();

% Center and radius of circular contour
c = 3i;
r = 2;

% Number of quadrature points
N = 500;

theta = linspace(0,2*pi,N+1);
theta = theta(1:end-1);
unit_circle = exp(1i*theta);

% Evenly-spaced quadrature points on circular contour
z = c + r*unit_circle;
% Quadrature weights
w = 2i*pi*r/N*unit_circle;

% Number of Hankel moments
k = 2;
% Number of probing directions
p = 7;

e = basic_solver(T,z,w,p,k)

