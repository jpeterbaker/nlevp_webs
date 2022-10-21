function [T,TV,gamma,nodes,edges] = regweb(spokes,rings,d)
%function [T,TV,gamma,nodes,edges] = regweb(spokes,rings,d)
%
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem
% The network is a "spider web" with regularly spaced spokes and rings
%
% For this problem
%    Nodes        nv = spokes*(rings+1)+1
%    Strings      ne = spokes*(2*rings+1)
%
% INPUTS
%
% spokes is the number of radial spokes in the web
% 
% rings is the number of concentric polygonal rings in the web
%
% d is the number of spatial dimensions for the web to occupy (2 or 3)
%     (default d=2)
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

if nargin<3 || isempty(d)
    d = 2;
end

[nodes,edges] = regweb_graph(spokes,rings,d);

[T,TV,gamma] = general_web(nodes,edges);

