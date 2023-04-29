% Create mesh
clear;
filename = strcat('data/laplace/problemToby/raw_meshes_square.mat');
load(filename);
elementType = 'P3gl';

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

filename = strcat('data/laplace/problemToby/meshSquare',elementType,'.mat');
save(filename,'mesh','dof');
%clear