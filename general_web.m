function [T,TV,gamma] = general_web(nodes,edges,rho,k,s,c)
%function [T,TV,gamma] = general_web(nodes,edges,rho,k,s,c)
%
% Represent a vibrating network of elastic strings as a nonlinear eigenvalue problem
%
% INPUTS
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
% rho is an ne-vector of string densities
%     rho(i) is the linear density of string i (when stretched in the network rest state)
%     If all edge densities are the same, rho may be a scalar
%     (default 1)
%
% k is an ne-vector of Hooke's constants for the strings
%     If all Hooke constants are the same, k may be a scalar
%     (default 1)
%
% s(i)>1 is the stretch factor of edge i (when web is at rest)
%     (default 2)
%
% c(i,1) is the viscous damping coefficient for string i in longitudinal direction
% c(i,2) is the viscous damping coefficient for string i in transverse directions
%     If c is a scalar, it is interpreted as the transverse damping for all strings
%     and longitudinal damping is zero.
%     (default 0)
%
% OUTPUTS
% 
% T is a function handle that accepts a scalar and returns a square matrix
% of dimension (2*d*ne) x (2*d*ne)
%     If T(lambda) is singular, then lambda/(2*pi) is a mode of the web
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
% gamma is an ne x 2 array of eigenvalues of the wave operator
%     gamma(i,1) is the eigenvalue for the longitudinal component of string i
%     gamma(i,2) for is the eigenvalue for all transverse directions of string i

[ne,two] = size(edges);
if two~=2
    error('edges must be ne x 2');
end
d = size(nodes,2);

if d>3
    warning('nodes should usually be (nv x 2) or (nv x 3)')
end

%--------------------------------%
% Combine fixed ends into node 1 %
%--------------------------------%
[node_degree,node_unique] = groupcounts(edges(:));
fixed = node_degree==1;
nv = numel(node_unique);

nf = sum(fixed);

u = zeros(nv,1);
u( fixed) = 1;
u(~fixed) = 2:nv-nf+1;

new_edges = zeros(ne,2);
for i=1:nv
    new_edges(edges==node_unique(i)) = u(i);
end

% Fill a sparse matrix all at once with lists of indices and entries
B = sparse([new_edges(:,1);new_edges(:,2)],[1:ne,1:ne]',[-ones(ne,1);ones(ne,1)]);

% From now on, nv is the number of nodes AFTER combining fixed ends
nv = nv - nf + 1;

%------------------------%
% Lengths and directions %
%------------------------%
L = zeros(ne,1);
vs = zeros(d,ne);
for i=1:ne
	tail = nodes(edges(i,1),:);
	head = nodes(edges(i,2),:);
	dx = head-tail;
	L(i) = norm(dx);
    if L(i) == 0
        error('Zero-length string not allowed')
    end
	vs(:,i) = dx/L(i);
end

if nargin < 4 || isempty(rho)
    rho = 1;
end
if nargin < 5 || isempty(k)
    k = 1;
end
if nargin < 6 || isempty(s)
    s = 2;
end
if nargin < 7 || isempty(c)
    c = 0;
end

if numel(rho) == 1
    rho = rho*ones(ne,1);
end
if numel(k) == 1
    k = k*ones(ne,1);
end
if numel(s) == 1
    s = s*ones(ne,1);
end

if numel(c) == 1
    c = [zeros(ne,1) c*ones(ne,1)];
end

if ~all(size(c) == [ne,2])
    error('damping c should be scalar or ne x 2');
end

L = L(:);
rho = rho(:);
k = k(:);
s = s(:);

TV = zeros(d,d,ne);
for i=1:ne
    % First column is string direction
    % Other columns are any orthonormal basis for the rest of R^d
    [TV(:,:,i),~] = svd(vs(:,i));
end

gamma = [ k.*s , k.*(s-1) ];

% build_mat was written for undamped problems and real frequencies
% Multiplying inputs by -1i changes the problem to the more standard convention
% where undamped eigenvalues are on the imaginary axis
% and damped eigenvalues are in the complex left half-plane.
S = @build_mat;
T = @(lambda) S(-1i*lambda);

function M = build_mat(omega)
g = bsxfun(@minus,omega^2*rho,1i*omega*c);
G = bsqrt(g);
pq = bsxfun(@times,L,G) ./ sqrt(gamma);
p = sin(pq);
q = cos(pq);

M = zeros(2*d*ne);
% First unused row
row0 = 1;

% Say frame has F string connections
% First d*F rows are Dirichlet conditions on the frame
% Next d*(2*ne-F-nv+1) rows are node continuity conditions
% Last d*(nv-1) rows are force balancing conditions

%------------------%
% Frame conditions %
%------------------%
Fi = find(B(1,:));
F = numel(Fi);
for i = Fi
    % String i is connected to the frame
    % so X and Y are  0 on this end
    % w=-1 if tail is connected to the frame, w=1 for head
    w = B(1,i);
    % First column associated with a coefficient for string i
    col0 = coeff_index(1,1,i);
    if w == -1
        % x=0 on the tail end, so cosine coefficients are zero
        % put 1 in matrix corresponding to all cosine coefficients
        for j=1:d
            M(row0+j-1 , col0+2*j-1) = 1;
        end
    else
        % First dimension
        M(row0,col0               ) = p(i,1);
        M(row0,col0+1             ) = q(i,1);
        % Other dimensions
        for j=2:d
            M(row0+j-1,col0+2*j-2:col0+2*j-1) = [p(i,2),q(i,2)];
        end
    end
    row0 = row0 + d;
end

%----------------------------%
% Node continuity conditions %
%----------------------------%
% for each non-frame node u
%     for all but the first connected edge
%         enter continuity condition: this edge and first edge agree
for u=2:nv
    % Edge-neighborhood of node u
    Nu = find(B(u,:));
    % First neighbor of u
    i1 = Nu(1);
    % Figure out position of string i1 at u
    w = B(u,i1);
    % yi1 will have the coefficients that would go in the matrix
    % if the string were horizontal.
    % TV(:,:,i1) has the correct transformation to standard coordinates
    % In other words, Y_i1(u) = yi1 * [A_i1 ; B_i1 ; ...]
    % and X_i1(u) = TV(:,:,i1)*Y_i1;
    yi1 = zeros(d,2*d);
    if w == -1
        % Node is at the tail (x=0), so sin terms are zero and cos terms are 1
        % If d=2, here's the result:
        % yi1 = [ 0 , 1 , 0 , 0 ; 0 , 0 , 0 , 1 ];
        for j=1:d
            yi1(j,2*j) = 1;
        end
    else
        % Node is at the head (x=L), so terms are already computed in p and q
        % First eig used for first direction, second eig for all others
        % If d=2, here's the result:
        % yi1 = [ p(i1,1) , q(i1,1)] , 1 , 0 , 0 ; 0 , 0 , p(i1,2) , q(i1,2) ];
        yi1(1,1:2) = [p(i1,1),q(i1,1)];
        for j=2:d
            yi1(j,2*j-1:2*j) = [p(i1,2),q(i1,2)];
        end
    end
    % -xi1 since it will be used like negxi1 + xi = 0
    negxi1 = -TV(:,:,i1)*yi1;
    col0i1 = coeff_index(1,1,i1);

    % Construction of yi is just like for yi1
    % except that, since it is reused, 0 must be written in w=-1 case
    yi = zeros(d,2*d);
    for i=Nu(2:end)
        w = B(u,i);
        if w == -1
            for j=1:d
                yi(j,2*j-1:2*j) = [0,1];
            end
        else
            yi(1,1:2) = [p(i,1),q(i,1)];
            for j=2:d
                yi(j,2*j-1:2*j) = [p(i,2),q(i,2)];
            end
        end
        xi = TV(:,:,i)*yi;

        col0i = coeff_index(1,1,i);
        M(row0:row0+d-1,col0i1:col0i1+2*d-1) = negxi1;
        M(row0:row0+d-1,col0i :col0i +2*d-1) = xi;
        row0 = row0 + d;
    end
end

%-------------------------------%
% Node force balance conditions %
%-------------------------------%
% Node 1 is the frame, but forces balance at all other nodes
% for each non-frame node u
%     for all connected edges
%         enter part of force condition: sum of forces is zero
% This is similar to the continuity loop
% but, for each node
%     only d rows are used (representing sum of all forces in each of d dimensions)
%     sin and cos (p and q) terms are modified (due to derivative)
%     a factor of sqrt( ( omega^2*rho + omega*c )*gamma ) appears
%         (due to derivative and multiplication by gamma matrix)
%
% The forces are H*dX but it's more convenient to write that as
%    TV*gamma*dY
% The inner loop calculates matrix entries corresponding to gamma*dY (Lyip)
% and then pre-multiplies by TV
for u=2:nv
    % Edge-neighborhood of node u
    Nu = find(B(u,:));

    % Y_i prime (derivative)
    Lyip = zeros(d,2*d);
    for i=Nu(1:end)
        w = B(u,i);
        % This term will be used repeatedly
        if w == -1
            Lyip(1,1:2) = [ sqrt(G(i)*gamma(i,1)) , 0 ];
            for j=2:d
                Lyip(j,2*j-1:2*j) = [ sqrt(G(i)*gamma(i,2)) , 0 ];
            end
        else
            Lyip(1,1:2) = sqrt(G(i)*gamma(i,1))*[q(i,1),-p(i,1)];
            for j=2:d
                Lyip(j,2*j-1:2*j) = sqrt(G(i)*gamma(i,2))*[q(i,2),-p(i,2)];
            end
        end
        Hxip = TV(:,:,i)*Lyip;

        col0i = coeff_index(1,1,i);
        % Factor of -w because derivative refers to different directions
        % depending on orientation
        M(row0:row0+d-1,col0i :col0i +2*d-1) = -w*Hxip;
    end
    row0 = row0 + d;
end
end % build_mat function

function ci = coeff_index(ab,j,i)
% Get the index of the coefficient
% ab = 1 for A, 2 for B
% d is dimension (in 1:d)
% i is string number
ci = ab+2*(j-1)+2*d*(i-1);
end % coeff_index function

end % main general_web function

