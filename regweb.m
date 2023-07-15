function [T,TV,gamma,nodes,edges] = regweb(spokes,rings,d)
%function [T,TV,gamma,nodes,edges] = regweb(spokes,rings,d)
%
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem.
% The network is a "spider web" with regularly spaced spokes and rings.
% The tension in each string is chosen to balance the forces at each free node.
%    Each ring strand has the same tension,
%    center spokes (connected to the hub) have tension 1 (stretch factor s=2),
%    and anchored spokes (on the perimeter) have tension 2 (s=3).
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

rho = 1;
k   = 1;

%--------------------------------------%
% Find force-balancing stretch factors %
%--------------------------------------%

% Inner angle between spoke and rings
theta = (0.5-1/spokes)*pi;
% Tension in each ring (scalar)
r_tension = sec(theta)/(2*rings);
% s_tension(i) is spoke tension between rings (i-1) and i
if rings==0
    s_tension = 1;
else
    s_tension = linspace(1,2,rings+1);
end
% Now replicate since there are (spokes) spoke-strings between consecutive rings
s_stretch = reshape(repmat(1+s_tension/k,spokes,1),[],1);

% Stretch factor s for all strings in the order provided by regweb_graph
% tension = k*(s-1)
s = [...
    s_stretch; ... % Spokes
    repmat(1+r_tension/k,rings*spokes,1)... % Rings
];

[T,TV,gamma] = general_web(nodes,edges,rho,k,s);

