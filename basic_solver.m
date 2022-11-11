function ew = basic_solver(T,z,w,p,k,tol)
%function ew = basic_solver(T,z,w,p,k,tol)
%
% Basic contour integration solver based on Beyn's algorithm
%
% 
% INPUTS
%
% T is a function handle that accepts a complex scalar and returns a square matrix.
%
% z is a complex vector specifying points on the integration contour.
%
% w is a real vector of quadrature weights corresponding to z.
%
% p is the number of probing directions to use.
%
% k is the number of moments to include in the Hankel matrices.
%     p*k should be no less than the number of eigenvalues in the contour.
%     More moments tend to be needed if eigenvalues in the contour have linearly
%     dependent eigenvectors.
%     (default k = 1)
%
% tol is the singular value cutoff for estimating the number of eigenvalues.
%     (default tol = 1e-6)
%
%
% OUTPUTS
%
% ew is a vector of approximate eigenvalues of T that lie in the contour defined by z.
%    ew is likely to be inaccurate if
%       > p or k is too small
%       > the contour points z and weights w do not correspond to a good quadrature method
%       > the contour points z are too coarse (with too few points)
%       > the contour is too large (with points that are far away from each other)
%       > the contour passes very close to an eigenvalue (inside or outside)
%

% by Jonathan Baker
% jonpbak@vt.edu

if nargin < 5
	k=1;
end
if nargin < 6
	tol = 1e-6;
end

% Size of the NLEVP
n = size(T(0),1);

N = numel(z);

% Set the random seed for reproducibility
rng(4);

L = randn(n,p);
R = randn(n,p);

A = zeros(p,p,2*k);

for j=1:N
    LTRj  = L'/T(z(j))*R;
    pA = zeros(p,p,2*k);
    for i=1:2*k
        % Add the contribution of T(j) to every A matrix
        Pij = w(j)*z(j)^(i-1)*LTRj;
        pA(:,:,i) = Pij;
    end
    A = A+pA;
end
% Because of the division that takes place when the B matrix is calculated,
% scaling A doesn't make a difference to the results
% except with respect to the singular value cutoff
A = A/(2i*pi);

H  = zeros(k*p);
Hs = zeros(k*p);

for i=1:2*k
    %%%%%%%%%%
    % Load H %
    %%%%%%%%%%
    % Number of matrices on the anti-diagonal of H containing Ai
    diaglen = k-abs(k-i);
    % Coordinates of top left corner of bottom left block
    % being filled in this iteration
    yhi = min(i-1,k-1)*p+1;
    xlo = max(0  ,i-k)*p+1;
    for j=1:diaglen
        inds = [yhi-(j-1)*p , yhi-(j-2)*p-1 , xlo+(j-1)*p , xlo+j*p-1];
        H(inds(1):inds(2),inds(3):inds(4)) = A(:,:,i);
    end
    %%%%%%%%%%%
    % Load Hs %
    %%%%%%%%%%%
    % Number of matrices on the anti-diagonal of Hs containing Ai
    diaglen = k-abs(k-i+1);
    % Indices of top left corner of bottom left block
    % being filled in this iteration
    yhi = min(i-2,k-1)*p+1;
    xlo = max(0,i-1-k)*p+1;
    for j=1:diaglen
        inds = [yhi-(j-1)*p , yhi-(j-2)*p-1 , xlo+(j-1)*p , xlo+j*p-1];
        Hs(inds(1):inds(2),inds(3):inds(4)) = A(:,:,i);
    end
end

[V0,S0,W0] = svd(H,0);
S0 = diag(S0);

% Estimate number of eigenvalues by singular values
ne = find(S0/S0(1)>tol,1,'last');

V0 = V0(:,1:ne);
S0 = S0(1:ne);
W0 = W0(:,1:ne);

%B = V0'*A1*W0/S0;
B = bsxfun(@rdivide,V0'*Hs*W0,S0.');
ew = eig(B);

