function [T,TV,gamma,nodes,edges,ew] = spider1(d)
%function [T,TV,gamma,nodes,edges] = spider1(d)
%
% Deinopis web model
% 
% Represent modal vibrations of a network of elastic strings as a nonlinear eigenvalue problem
% The network is a tracing of the web of a Deinopis spider
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

load('deinopis.mat')

% If a third dimension is requested, put web in z=0 plane
if d == 3
    nodes(1,3) = 0;
end

[T,TV,gamma] = general_web(nodes,edges);

ew = 1i*[
    0.00430521
    0.00500223
    0.00607796
    0.00684835
    0.00775994
    0.00808546
    0.00920314
    0.00967117
    0.01120355
    0.01171892
    0.01295983
    0.01319632
    0.01341148
    0.01411616
    0.01527657
    0.01625177
    0.01688942
    0.01732988
    0.01786737
    0.01794647
    0.01856216
    0.01923194
    0.01984739
    0.02019651];

