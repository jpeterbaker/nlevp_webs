
% This script demonstrates how to produce a string NLEVP and visualize a mode shape

%-------------------------------------%
% Index of the mode to display (1-24) %
% Please experiment with this value   %
%-------------------------------------%
j = 7;

% Produce the matrix function (T),
%     network description (nodes, edges)
%     supplemental information for plotting (TV, gamma),
%     and the first 24 eigenvalues (ew)
[T,TV,gamma,nodes,edges,ew] = tritare(pi/2);

% Select one eigenvalue for study
lambda = ew(j);

% Get an associated eigenvector v
[~,~,V] = svd(T(lambda));
v = V(:,end);

% Get string positions during modal vibration
R = mode_curves(lambda,v,TV,gamma,nodes,edges,1,0,0.3,60);

clf
hold on

% Plot at-rest network in grey
plot([nodes(edges(:,1),1)' ; nodes(edges(:,2),1)'],[nodes(edges(:,1),2)' ; nodes(edges(:,2),2)'],'-','color',[0.7,0.7,0.7]);

% Plot mode shape in black with dots for visibility of longitudinal motion
for i=1:3
    plot(R(1,:,i),R(2,:,i),'k.-',...
        'markersize',16,...
        'markerindices',1:3:60 ...
    );

end
title(sprintf('Mode shape of $\\omega_{%d}\\approx%0.5f$',j,imag(lambda)),'interpreter','latex');

axis off
axis equal

