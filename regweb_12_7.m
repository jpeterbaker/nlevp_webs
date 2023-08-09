function [T,TV,gamma,nodes,edges,ew] = regweb_12_7(d)
%function [T,TV,gamma,nodes,edges,ew] = regweb_12_7(d)
%
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem
% The network is a "spider web" with 12 spokes and 7 rings
%
% For this problem
%    Nodes        nv = 97
%    Strings      ne = 180
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
[T,TV,gamma,nodes,edges] = regweb(12,7,d);

ew = 1i*[
    1.554718968135
    1.554718968135
    2.392757029983
    2.733524641834
    2.733524641834
    2.972047754465
    3.198512077458
    3.198512077458
    3.286656934694
    3.286656934694
    3.494873518374
    3.494873518374
    3.526354924854
    3.526354924854
    3.553608948641
    3.553608948641
    3.571196078490
    3.597944965839
    3.597944965839
    3.624703655060
    3.658276748550
    3.658276748550
    3.741538406894
    3.741538406894];

