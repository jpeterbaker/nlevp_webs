function [B,new_edges,vs,L] = incidence_vectors(nodes,edges)
%function [B,new_edges,vs,L] = incidence_vectors(nodes,edges)
%
% Build the "incidence matrix" of the network
% and calculate the direction and length of each string.
% Nodes with degree 1 are presumed to be fixed anchors.
% In the incidence matrix, the fixed string ends are consolidated into a single node.
%
% INPUTS
%
% nodes is an nv x n matrix of node coordinates
%     nodes(i,:) is the spatial location of node i
%
% edges is an ne x 2 matrix of node indices
%     for each i, edge i is directed from node edge(i,1) to edge(i,2)
%     maximum value in edges is nv
%
% OUTPUTS
%
% B is an (nv-nf) x ne matrix
%     nf is the number of fixed ends (nodes with degree 1)
%     B(i,j) is 1 if node i is the head of edge j
%              -1 if node i is the tail of edge j
%               0 otherwise
%
% new_edges is an (nv-nf) x 2 matrix of node indices
%     It is comparable to the original edges matrix except that
%     the anchored ends of strings have been combined into the node with index 1
%
% vs is an n x ne matrix of edge orientations
%     vs(:,i) is a unit vector in the direction from the tail to the head of edge i
%
% L is an ne-vector of edge lengths
%     L(i) is the length of edge i

[ne,two] = size(edges);
if two~=2
    error('edges must be ne x 2');
end
n = size(nodes,2);

% Find fixed ends and combine their nodes into node 1
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

L = zeros(ne,1);
vs = zeros(n,ne);
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

