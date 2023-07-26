function [nodes,edges] = regweb_graph(ns,nr,d)
%function [nodes,edges] = regweb_graph(ns,nr,d)
%
% Generate the node positions and edge connections in a
% highly regular prototypical spider web
%
% INPUTS
%
% ns is the number of radial spokes in the web
% 
% nr is the number of concentric polygonal rings in the web
%
% d is the number of spatial dimensions for the web to occupy (2 or 3)
%     (default  d=2)
%
% OUTPUTS
% 
% nodes is an nv x 2 matrix with the coordinates of nodes
%     nodes(i,:) is the location of node i
%     nv is the number of nodes: ns*(nr+1)+1
% 
% edges is an ne x 2 matrix with the indices of connected nodes
%     edges(i,:) contains the indices of the two nodes connected by string i
%     ne is the number of edges: ns*(2*nr+1)
%
if nargin<3 || isempty(d)
    d=2;
end
if d==1
    error('d must be >= 2')
end

% Angles at which to place spokes
theta = linspace(0,2*pi,ns+1);
theta = theta(1:end-1);

% Radii at which to place rings
rs = linspace(0,1,nr+2);
rs = rs(2:end-1);

nodes = [0;reshape( exp(1i*theta).'*[rs 1], [],1)];
nodes = [real(nodes) imag(nodes)];

nv = ns*(nr+1)+1;
ne = ns*(2*nr+1);

edges = zeros(ne,2);

% Edges of center node
edges(1:ns,:) = [ones(ns,1) (2:ns+1)'];
% Spokes
edges(ns+1:nv-1,1) = 2:nv-ns;
edges(ns+1:nv-1,2) = 2+ns:nv;
% Rings
edges(nv:2*nv-ns-2,1) = 2:nv-ns;
edges(nv:2*nv-ns-2,2) = (3:nv-ns+1)-ns*(mod(2:nv-ns,ns)==1);

nodes(1,d) = 0;

