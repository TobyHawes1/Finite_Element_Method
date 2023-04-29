%To run code:
%pathmd
%beginlabfem2d(@pathmd)
%beginproblem(@pathmd,'laplace/problemToby/');
%executeCode

clear;
close all;

elementType = ['P3_3'];

% Load mesh
filename = strcat('data/laplace/problemToby/meshSquare',elementType,'.mat');
% elementType = 'P3gl';
% filename = strcat('data/laplace/problemToby/meshSquare',elementType,'.mat');
load(filename);
clear elementType filename;

% Load settings and solve problem
settings = loadSettings;

[u,error] = laplace(mesh,dof,settings);