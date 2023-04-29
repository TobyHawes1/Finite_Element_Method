%load('data/laplace/problemToby/meshP1_1.mat');

% P1 ELEMENTS
% Plot mesh and highlight P1 nodes
%load('data/laplace/problemToby/meshP1_1.mat');
p = p1;
e = e1;
t = t1;
close all
figure(1);
pdeplot(p,t(1:3,:));
hold on;
plot(mesh.d1_v1_p(1,:),mesh.d1_v1_p(2,:),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
% Plot one basis function
node = randi(dof.d1_nv1);
u.d1_v1 = zeros(dof.d1_nv1,1);
u.d1_v1(node,1) = 1;
plotLaplace(mesh,'d1',u,0,0,1);
title('Global basis function for P1 polynomials')
%

% P2 ELEMENTS
% Plot mesh and highlight P2 nodes
clear mesh dof;
load('data/laplace/problemToby/meshSquareP2.mat');
figure(3);
pdeplot(p,t(1:3,:));
hold on;
plot(mesh.d1_v1_p(1,:),mesh.d1_v1_p(2,:),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
% Plot one basis function (still to do)
%

% P3 ELEMENTS
% Plot mesh and highlight P3 nodes
clear mesh dof;
load('data/laplace/problemToby/meshSquareP3.mat');
figure(4);
pdeplot(p,t(1:3,:));
hold on;
plot(mesh.d1_v1_p(1,:),mesh.d1_v1_p(2,:),'o','MarkerEdgeColor','b','MarkerFaceColor','b');
% Plot one basis function (still to do)
%