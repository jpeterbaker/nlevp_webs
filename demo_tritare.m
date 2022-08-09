
% Demonstrates how to produce and visualize a string NLEVP

% Produce the matrix function (T),
%     network description (nodes, edges)
%     supplemental information for plotting (TV, gamma),
%     and the first 24 positive eigenvalues (eigs)
[T,TV,gamma,nodes,edges,eigs] = tritare();

% Select one eigenvalue for study
omega = eigs(3);
T(omega)

% Get an associated eigenvector v
[~,~,V] = svd(T(omega));
v = V(:,end);

% Get string positions during modal vibration
R = mode_curves(omega,v,TV,gamma,nodes,edges,[],[],0.2);

fig=figure(1);
clf
hold on

% Plot at-rest network in grey
plot([nodes(edges(:,1),1)' ; nodes(edges(:,2),1)'],[nodes(edges(:,1),2)' ; nodes(edges(:,2),2)'],'-','color',[0.7,0.7,0.7]);

% Plot mode shape in black
for i=1:3
    plot(R(1,:,i),R(2,:,i),'k-');
end
title(['Mode shape of $\omega$=',sprintf('%g',omega)],'interpreter','latex');

axis off


