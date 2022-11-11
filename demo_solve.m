function ew = demo_solve(num)
%function ew = demo_solve(num)
%
% Demonstrates use of the basic integration solver, basic_solver.m,
% using a selection of the vibrating string problems in this package.
%
% The selected string network will be drawn in a figure along with
% some of its eigenvalues and the (predetermined) contour used to calculate them.
%
% The eigenvalues tend to be inaccurate if the contour is too large or passes
% too near an eigenvalue. Eigenvalues that lie near but outside the contour 
% are sometimes calculated unexpectedly by this method.
% The contours in this demo were selected to avoid these issues.
% Contours that are chosen with less care may not give such good results.
% 
% INPUT
%
% num is a problem selector
% 1) tritare (Y-shaped string)
% 2) regweb_5_4 (prototypical web with 5 spokes, 4 rings)
% 3) regweb_12_7 (prototypical web with 12 spokes, 7 rings)
% 4) spider1 (Deinopis web)
% 5) spider2 (Orb weaver web)
%
% If num is not provided as a parameter,
% user is prompted to enter it manually.
%
% OUTPUT
%
% ew is a vector of the calculated eigenvalues
% 

% by Jonathan Baker
% jonpbak@vt.edu

%-------------%
% Check input %
%-------------%
if nargin < 1
    newline = sprintf('\n');
    helpstr = join([
        "Enter a problem number to continue:",
        "1) tritare (Y-shaped string)",
        "2) regweb_5_4 (prototypical web with 5 spokes, 4 rings)",
        "3) regweb_12_7 (prototypical web with 12 spokes, 7 rings)",
        "4) spider1 (Deinopis web)",
        "5) spider2 (Orb weaver web)",
        ""],...
    newline...
    );
    num = input(helpstr,'s');
    num = str2double(num);
end

%----------------------------%
% Get selected NLEVP problem %
% and good predetermined     %
% solver parameters          %
%----------------------------%

% Number of quadrature points
N = 100;

if num == 1
    % The NLEVP function
    [T,~,~,nodes,edges] = tritare();

    % Establish extent of elliptical contour
    xlo = -0.5; xhi =  0.5;
    ylo =  1.0; yhi = 10.0;

    % Number of Hankel moments
    k = 3;
    % Number of probing directions
    p = 10;

    fprintf("\nSolving tritare problem\n")
    fprintf("This will probably take less than a second\n\n")
elseif num == 2
    [T,~,~,nodes,edges] = regweb_5_4();

    xlo = -0.1; xhi = 0.1;
    ylo =  1.0; yhi = 4.6;

    k = 1;
    p = 30;

    fprintf("\nSolving 5-spoke 4-ring web problem\n")
    fprintf("This will probably take less than a second\n\n")
elseif num == 3
    [T,~,~,nodes,edges] = regweb_12_7();

    xlo = -0.1; xhi = 0.1;
    ylo =  1.0; yhi = 5.0;

    k = 1;
    p = 20;

    fprintf("\nSolving 12-spoke 7-ring web problem\n")
    fprintf("This will probably take a few seconds\n\n")
elseif num == 4
    [T,~,~,nodes,edges] = spider1();

    xlo = -0.005; xhi = 0.005;
    ylo =  0.001; yhi = 0.0203;

    k = 1;
    p = 25;

    fprintf("\nSolving net-caster web problem\n")
    fprintf("This will probably take less than a second\n\n")
elseif num == 5
    [T,~,~,nodes,edges] = spider2();

    xlo = -0.005; xhi = 0.005;
    ylo =  0.001; yhi = 0.009;

    k = 1;
    p = 10;

    fprintf("\nSolving orb-weaver web problem\n")
    fprintf("This will probably take a few minutes\n\n")
else
    error("Input should be an integer 1-5");
end

%------------------------------------------------%
% Set up points on contour and quadrature points %
%------------------------------------------------%
theta = linspace(0,2*pi,N+1);
theta = theta(1:end-1);
a = (xhi-xlo)/2;
b = (yhi-ylo)/2;
c = (xhi+xlo)/2 + 1i*(yhi+ylo)/2;
% Points on contour
z  = c + a* cos(theta) + 1i*b*sin(theta);
% Quadrature weights
w = 2*pi/N*(a*-sin(theta) + 1i*b*cos(theta));

%------------------%
% Use NLEVP solver %
%------------------%
ew = basic_solver(T,z,w,p,k);

%----------------------------------%
% Plot the contour and eigenvalues %
%----------------------------------%
figure(10)
subplot(1,2,1)
hold on

plot(real([z,z(1)]),imag([z,z(1)]),'b-')
plot(real(ew),imag(ew),'ko')
axis equal
title('Eigenvalues and integration contour')

%------------------%
% Draw the strings %
%------------------%
subplot(1,2,2)
hold on

% Number of edges
ne = size(edges,1);

for i=1:ne
	ei = edges(i,:);
	plot(nodes(ei,1),nodes(ei,2),'k-');
end
axis equal
title('String network at rest')

