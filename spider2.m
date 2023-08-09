function [T,TV,gamma,nodes,edges,ew] = spider2(d)
%function [T,TV,gamma,nodes,edges,ew] = spider2(d)
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
% ew is a vector of the first 24 eigenvalues
%     ONLY AVAILABLE IF d IS 2
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

if d == 2
    ew = 1i*[
        0.02683803
        0.02872040
        0.06095035
        0.06973589
        0.07226500
        0.07949915
        0.09029340
        0.09669386
        0.10341368
        0.10521216
        0.11192470
        0.11679917
        0.11729910
        0.12723586
        0.13377063
        0.13767223
        0.13956348
        0.14561710
        0.14608491
        0.15123553
        0.15362060
        0.15819261
        0.16258235
        0.16368702];
end

