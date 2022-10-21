
% This script demonstrates how to produce and visualize a string NLEVP

%-------------------------------------%
% Index of the mode to display (1-24) %
% Please experiment with this value   %
%-------------------------------------%
j = 3;

% Produce the matrix function (T),
%     network description (nodes, edges)
%     supplemental information for plotting (TV, gamma),
%     and the first 24 eigenvalues (eigs)
[T,TV,gamma,nodes,edges,eigs] = tritare();

% Select one eigenvalue for study
lambda = eigs(j);

% Get an associated eigenvector v
[~,~,V] = svd(T(lambda));
v = V(:,end);

% Get string positions during modal vibration
R = mode_curves(lambda,v,TV,gamma,nodes,edges,[],[],0.2);

figure(1)
clf
hold on

% Plot at-rest network in grey
plot([nodes(edges(:,1),1)' ; nodes(edges(:,2),1)'],[nodes(edges(:,1),2)' ; nodes(edges(:,2),2)'],'-','color',[0.7,0.7,0.7]);

% Plot mode shape in black
for i=1:3
    plot(R(1,:,i),R(2,:,i),'k-');
end
title(['Mode shape of $\lambda$=',sprintf('%gi',lambda/1i)],'interpreter','latex');

axis off

