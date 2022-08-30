function R = mode_curves(lambda,v,TV,gamma,nodes,edges,rho,c,scale,N)
%function R = mode_curves(lambda,v,TV,gamma,nodes,edges,rho,c,scale,N)
%
% Get string vibration positions for a particular dynamic mode in d dimensions
%
% INPUTS
%
% lambda is an eigenvalue of the NLEVP
%
% v is an eigenvector associated with lambda
%     i.e. a singular vector of T(lambda)
%
% TV is the array of basis matrices produced by the NLEVP generation functions
%
% gamma is the matrix of wave operator eigenvalues produced by the NLEVP generation functions
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
% rho is a scalar or ne-vector of linear densities of the strings
%     default 1
%
% c is a scalar or ne x d matrix of damping coefficients
%     default 0
%
% scale is the mode displacement multiplier
%     default 1
%
% N is the number of points with which to discretize each string
%     default 30 
%
% OUTPUT
%
% R is an d x N x ne array
%     R(:,i,j) is is the position of the ith bit of string j
%

if nargin < 7 || isempty(rho)
    rho = 1;
end
if nargin < 8 || isempty(c)
    c = 0;
end
if nargin < 9 || isempty(scale)
    scale = 1;
end
if nargin < 10 || isempty(N)
    N = 30;
end

if max(abs(imag(v(:)))) > 1e-6
    error('v should be a real vector. Find a real basis for span of eigenvalue''s null space.')
end

v = real(v);

[nv,d] = size(nodes);
[ne,two] = size(edges);

if two ~= 2
    error('edges must be ne x 2');
end

if d>2 && size(gamma,2)<d
    % Add duplicate columns to gamma if they were omitted
    gamma = [gamma , repmat(gamma(:,end),1,d-size(gamma,2))];
end

if numel(rho) == 1
    rho = rho*ones(ne,1);
end
if numel(c) == 1
    c = c*ones(ne,1);
end

[~,~,vs,L] = incidence_vectors(nodes,edges);

% All of the tail locations
R0 = nodes(edges(:,1),:)';

R = zeros(d,N,ne);

% This term gets multiplied by eta and then becomes input to sin or cos
PQ = sqrt(bsxfun(@rdivide,-bsxfun(@plus,lambda^2*rho,lambda*c),gamma));
for i=1:ne
    eta = linspace(0,L(i),N);
    % pq for this string in all d dimensions
    pq = PQ(i,:)';
    pqeta = bsxfun(@times,pq,eta);
    % Two entries per dimension per string
    k = (i-1)*2*d;
    Y = ...
        bsxfun(@times,v(k+1:2:k+2*d),sin(pqeta)) + ...
        bsxfun(@times,v(k+2:2:k+2*d),cos(pqeta));
    X = TV(:,:,i)*Y;
    % Static equilibrium (rest) position of the stretched string
    rest = bsxfun(@plus,R0(:,i), bsxfun(@times,eta,vs(:,i)) );
    R(:,:,i) = rest + scale*X;
end

