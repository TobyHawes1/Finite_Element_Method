% Create mesh
clear;
filename = strcat('data/laplace/problemToby/Sparse_Mesh_Nodes_Example.mat');
load(filename);
elementType = 'P3gl';
whichMesh = 3;

switch whichMesh
    case 1
        p = p1;
        e = e1;
        t = t1;
    case 2
        p = p2;
        e = e2;
        t = t2;
    case 3
        p = p3;
        e = e3;
        t = t3;
    case 4
        p = p4;
        e = e4;
        t = t4;
end

switch elementType
    case 'P1'
        [mesh] = meshP(p,e,t,1);
    case 'P2'
        [mesh] = meshP(p,e,t,2);
    case 'P3'
        [mesh] = meshP(p,e,t,3);
    case 'P3gl'
        [mesh] = meshP(p,e,t,31);
end

variableName = 'v1';
variableComponents = 1;
domainName = 'd1';
infoProblemMesh(1,:) = {mesh,variableName,variableComponents,elementType,domainName};
[mesh,dof] = prepareMesh (infoProblemMesh,[],[]);

filename = strcat('/Users/tobyhawes/Documents/Dissertation/labfem2d-3/data/laplace/problemToby/meshTest',elementType,'_',int2str(whichMesh),'.mat');
save(filename,'mesh','dof');
%clear