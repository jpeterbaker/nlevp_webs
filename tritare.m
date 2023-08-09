function [T,TV,gamma,nodes,edges,ew] = tritare(theta)
%function [T,TV,gamma,nodes,edges,ew] = tritare(theta)
%
%    If theta is one of 2*pi/3 or pi/2, then pre-computed eigenvalues ew
%    are available with the call signature
% 
%function [T,TV,gamma,nodes,edges,ew] = tritare(theta)
%
%
% Produces an NLEVP corresponding to planar modal vibrations
% of a Y-shaped "three-string" or "tritare"
%
% For each string, Hooke's constant (k) and linear density (rho) are 1.
% Tension in each "arm" has magnitude 1.
% Tension in the "stem" is selected to produce static equilibrium
% when the arms form the angle theta.
%
% INPUT
%     theta: angle between the "arm" strings.
%         theta should be in the interval [0,pi).
%             (default value 2*pi/3)
%         The "stem" string lies opposite the bisector of this angle
%         so that it lies on a line of symmetry.
%
% OUTPUTS
% 
% T is a function handle that accepts a scalar and returns a square matrix
% of dimension (2*d*ne) x (2*d*ne)
%     If T(omega) is singular, then omega/(2*pi) is a natural frequency of the web
%     The corresponding singular vector gives coefficients of mode shape
%     The mode shape is described as a function X(x) for each string
%     Y(x) is X(x) transformed to a coordinate system parallel to the string
%     Each Y(x) has the form A*sin(a*x)+B*cos(b*x) in each spatial dimension
%
%     The coefficient order is
%     A for string 1 dimension 1
%     B for string 1 dimension 1
%     A for string 1 dimension 2
%     ...
%     A for string 1 dimension d
%     B for string 1 dimension d
%     A for string 2 dimension 1
%     ...
%     B for string nv dimension d
%
%     a,b for each string and dimension are determined by physical parameters
%
% TV is an d x d x ne array of basis vectors
%     This is needed for understanding mode shapes
%     TV(:,i,j) is the direction of "dimension i" for string j
%
% gamma is an ne x d array of eigenvalues of the wave operator
%     gamma(i,1) is the eigenvalue for the longitudinal component of string i
%     gamma(i,j) for j=2..d is the eigenvalue for all transverse directions of string i
%          The same value is repeated d-1 times
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
%     ONLY AVAILABLE IF theta IS 2*pi/3 OR pi/2
%

if nargin < 1
    theta = 2*pi/3;
end

% Magnitude of tension in arms (strings 2 and 3)
T2 = 1;
% Magnitude of tension in stem (string 1)
T1 = 2*cos(theta/2);
% Stretched string lengths (each has relaxed length 1 and stretch factor T+1)
L1 = T1+1;
L2 = T2+1;

nodes = [  0              ,   0
         -L1              ,   0
          L2*cos(theta/2) ,  L2*sin(theta/2)
          L2*cos(theta/2) , -L2*sin(theta/2)];

edges = [2,1
         3,1
         4,1];

[T,TV,gamma] = general_web(nodes,edges,1,1,[L1,L2,L2],0);

if nargout >= 5
    if theta == pi/2
        ew = 1i*[
        0.883875509098001370945688693572
        0.969600027367207098059069481547
        1.559940164853354585952582011179
        1.687255902893433860920888107107
        1.889068023025079978677377906620
        2.124080281919104993088069907433
        2.683784414143711570757257426799
        2.856805167333002858142965926867
        3.120353761893815338169926503391
        3.449638737120018841876201197701
        3.713034093361836938217116163254
        4.229510738816253833968713911968
        4.550207155201949550853215272526
        4.633725274842637194736881793762
        4.684483850268532793755600842737
        5.344027830827305043568241527448
        5.484411832262401894731713119148
        6.196616841084675471324621360538
        6.236185824688174223474712494009
        6.458002404821112841837438349256
        6.463341451576820650948480144456
        7.264515230048113874777574654440
        7.307505374754784315960662732216
        7.800149654953724042768337409983];
    elseif theta == 2*pi/3
        ew = 1i*[
        0.946647610817252850956049097461
        0.946647610817252850956049097461
        1.570796326794896619231321691640
        1.796490074799357190552357521150
        1.796490074799357190552357521150
        2.221441469079183123507940495030
        2.810809735860061632799891901174
        2.810809735860061632799891901174
        3.141592653589793238462643383280
        3.637882729207750660693570426221
        3.637882729207750660693570426221
        4.442882938158366247015880990061
        4.621912020675813962567144032789
        4.621912020675813962567144032789
        4.712388980384689857693965074919
        5.525688165978060107380202340099
        5.525688165978060107380202340099
        6.283185307179586476925286766559
        6.412038049334510793827550405036
        6.412038049334510793827550405036
        6.664324407237549370523821485091
        7.408212465121425903065818408268
        7.408212465121425903065818408268
        7.853981633974483096156608458199];
    else
        error('Precomputed eigenvalues only available for theta = pi/2 and 2*pi/3');
    end
end
