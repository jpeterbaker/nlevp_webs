function [T,TV,gamma] = regweb_4_3(d)
%function [T,TV,gamma] = regweb_4_3(d)
%
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem
% The network is a "spider web" with 4 spokes and 3 rings
%
% For this problem
%    Nodes        nv = 17
%    Strings      ne = 28
%
% INPUTS
%
% d is the number of spatial dimensions for the web to occupy (2 or 3)
%     default 2
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
if nargin<1 || isempty(d)
    d=2;
end

[T,TV,gamma,nodes,edges] = regweb(4,3,d);

