function [nodes,edges] = regweb_graph(spokes,rings,d)
%function [nodes,edges] = regweb_graph(spokes,rings,d)
%
% Generate the node positions and edge connections in a
% highly regular prototypical spider web
%
% INPUTS
%
% spokes is the number of radial spokes in the web
% 
% rings is the number of concentric polygonal rings in the web
%
% d is the number of spatial dimensions for the web to occupy (2 or 3)
%     default 2
%
% OUTPUTS
% 
% nodes is an nv x 2 matrix with the coordinates of nodes
%     nodes(i,:) is the location of node i
%     nv is the number of nodes: spokes*(rings+1)+1
% 
% edges is an ne x 2 matrix with the indices of connected nodes
%     edges(i,:) contains the indices of the two nodes connected by string i
%     ne is the number of edges: spokes*(2*rings+1)
%
if nargin<3 || isempty(d)
    d=2;
end
if d==1
    error('d must be >= 2')
end

% Angles at which to place spokes
theta = linspace(0,2*pi,spokes+1);
theta = theta(1:end-1);

% Radii at which to place rings
rs = linspace(0,1,rings+2);
rs = rs(2:end-1);

nodes = [0;reshape( exp(1i*theta).'*[rs 1], [],1)];
nodes = [real(nodes) imag(nodes)];

nv = spokes*(rings+1)+1;
ne = spokes*(2*rings+1);

edges = zeros(ne,2);

% Edges of center node
edges(1:spokes,:) = [ones(spokes,1) (2:spokes+1)'];
% Spokes
edges(spokes+1:nv-1,1) = 2:nv-spokes;
edges(spokes+1:nv-1,2) = 2+spokes:nv;
% Rings
edges(nv:2*nv-spokes-2,1) = 2:nv-spokes;
edges(nv:2*nv-spokes-2,2) = (3:nv-spokes+1)-spokes*(mod(2:nv-spokes,spokes)==1);

nodes(1,d) = 0;

