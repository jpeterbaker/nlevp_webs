function ew = demo_solve(num)
%function ew = demo_solve(num)
%
% Demonstrates use of the basic integration solver, basic_solver.m,
% using a selection of the vibrating string problems in this package.
%
% The selected string network will be drawn in a new figure along with
% some of its eigenvalues and the (predetermined) contour used to calculate them.
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
if nargin < 1 || numel(num) ~= 1
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
    num = str2num(num);
end

if num == 1
    % The NLEVP function
    T = tritare();

    % Center and radius of circular contour
    c = 3i;
    r = 2;

    % Number of quadrature points
    N = 500;

    % Set up evenly-spaced quadrature points on circular contour
    theta = linspace(0,2*pi,N+1);
    theta = theta(1:end-1);
    unit_circle = exp(1i*theta);
    z = c + r*unit_circle;
    % Quadrature weights
    w = 2i*pi*r/N*unit_circle;

    % Number of Hankel moments
    k = 2;
    % Number of probing directions
    p = 7;

    % Use NLEVP solver
    e = basic_solver(T,z,w,p,k)

    % Plot the contour and eigenvalues
    figure()
    hold on
    plot(real(z),imag(z),'k-')
    plot(real(e),imag(e),'ko')
    axis equal
elseif num == 2
elseif num == 3
elseif num == 4
elseif num == 5
else
    error("Input should be an integer 1-5");
end

