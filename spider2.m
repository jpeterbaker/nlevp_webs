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
% See general_web.m for an explanation of T,TV,gamma
%

if nargin<3 || isempty(d)
    d = 2;
end

load('orb.mat')

% If a third dimension is requested, put web in z=0 plane
if d == 3
    nodes(1,3) = 0;
end

[T,TV,gamma] = general_web(nodes,edges);

