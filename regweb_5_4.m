function [T,TV,gamma,nodes,edges,ew] = regweb_5_4(d)
%function [T,TV,gamma,nodes,edges,ew] = regweb_5_4(d)
%
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem
% The network is a "spider web" with 5 spokes and 4 rings
%
% For this problem
%    Nodes        nv = 26
%    Strings      ne = 45
%
% INPUTS
%
% d is the number of spatial dimensions for the web to occupy (2 or 3)
%     (default  d=2)
%
% OUTPUTS
% 
% See general_web.m for an explanation of T,TV,gamma
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
% ew is a vector of the first 24 eigenvalues
%
if nargin<1 || isempty(d)
    d=2;
end

[T,TV,gamma,nodes,edges] = regweb(5,4,d);

ew = 1i*[
    1.312358145733
    1.312358145733
    1.495788032320
    1.504476353083
    1.504476353083
    1.568414224554
    1.568414224554
    1.943739604204
    1.953508261143
    1.953508261143
    2.104406671255
    2.104406671255
    2.156096894499
    2.685496149881
    2.685496149881
    2.831708342636
    2.890933407419
    2.890933407419
    2.947006361942
    2.947006361942
    3.028297450748
    3.084589981280
    3.084589981280
    3.217829935846];

