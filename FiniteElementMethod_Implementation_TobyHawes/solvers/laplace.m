function [u,error,varargout] = laplace (mesh,dof,settings,varargin)

% laplace.m solves the Poisson problem in 2D using finite elements method.
% 
% INPUT:
%
%    mesh      structure for either P1 or P2 elements
%              created using .../Meshes/SetMeshScalar.m
%    dof       structure with the degrees of freedom
%              created using .../Meshes/SetMeshScalar.m
%    settings  structure created by LoadSettings.m
%              .method   {Direct}/PCG/PCG0
%              .NI       1/{0}
%              .plot     {1}/0
%              .error    1/{0}
%
% DATA OF THE PROBLEM:
%    To specify the data of the problem, modify the file LaplaceData.
%
% OUTPUT:
%
%    U structure containing the solution
%
%    ERROR is an optional structure containing the errors with respect to
%          the exact solution (if known) in norm L2 and H1.
%
% [u,error] = Laplace(mesh,dof,settings);
%

domain = 'd1';
variable = 'v1';
time = 0;

% Assemble the matrices
% ( Possible modifications due to Robin boundary conditions are made in ImposeLaplaceBC )
switch settings.NI
    case 0
        whichMatrix = [1,0,1,0];
    case 1
        whichMatrix = [-1,0,1,0];
end
[A,rhs,matrixData] = laplaceMatrices(mesh,dof,domain,variable,time,@laplaceData,{},whichMatrix);

% Impose boundary conditions
[dof] = selectBoundaryLaplace(dof,mesh,laplaceData(20),domain,variable);
[A,rhs,dirichletValues] = imposeLaplaceBC(A,matrixData,rhs,mesh,dof,domain,variable,time,@laplaceData,{});

% Modification for periodic boundary conditions (if used)
boundaryData = laplaceData(20);
periodicEdges = strcat(domain,'_PeriodicEdges');
if isfield(boundaryData,periodicEdges)
    [A,rhs,dof] = imposeLaplacePeriodicBC (A,matrixData,rhs,dof,domain,variable);
end
%

[A,rhs] = setLaplaceSystem(A,matrixData,rhs,dof,domain,variable);

% Solver
methodFlag = 0;
[matrixData] = allocateLaplaceMatrix(matrixData,settings,methodFlag,domain,[]);
domainVariable = strcat(domain,'_',variable);
[u] = solveLaplace(A,matrixData,rhs.(domainVariable),settings,methodFlag);
if nargout==3
    varargout(1) = {A.(matrixData.C)};
end
clear A rhs;

% Set boundary conditions
[u] = setLaplaceBC(u,dirichletValues,dof,domain,variable);
% Set periodic boundary conditions (if used)
if isfield(boundaryData,periodicEdges)
    [u] = setLaplacePeriodicBC (u,dof,domain,variable);
else
    clear boundaryData periodicEdges;
end
%

% Plot of the solution
if settings.plot
    plotLaplace (mesh,domain,u,settings.plotExact,time,[]);
end

% Compute error estimates (if required)
if settings.error
    error = laplaceError(u,mesh,domain,variable,@laplaceData,time);
else
    error = [];
end

end