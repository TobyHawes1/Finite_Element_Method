function [u,error] = laplaceTime (mesh,dof,t0,tf,dt,settings,varargin)

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
domainVariable = strcat(domain,'_',variable);
numberOfNodes = strcat(domain,'_n',variable);
nodes = strcat(domainVariable,'_p');
bcType = strcat(domain,'_',variable,'_bctype');
internalIndices = strcat(domainVariable,'_1_internal');
numberOfInternalNodes = strcat(domainVariable,'_1_ninternal');

% Initialisation
t = (t0:dt:tf)';
N = size(t,1);
u.(domainVariable) = zeros(dof.(numberOfNodes),N);
x = mesh.(nodes)(1,:);
y = mesh.(nodes)(2,:);
dataFlagInitialCondition = 40;
u.(domainVariable)(:,1) = laplaceData(dataFlagInitialCondition,x,y,t0)';

% Assemble the stiffness and mass matrices
switch settings.NI
    case 0
        whichMatrix = [ 1,1,0,0];
    case 1
        whichMatrix = [-1,1,0,0];
end
[A,rhs,matrixData] = laplaceMatrices(mesh,dof,domain,variable,[],@laplaceData,{},whichMatrix);

% Check boundary conditions
[dof] = selectBoundaryLaplace(dof,mesh,laplaceData(20),domain,variable);

% Prepare additional matrices depending on imposed boundary conditions
% Handle possible Robin boundary conditions
if dof.(bcType)(3)
    % Modify the stiffness matrix
    matrixData.boundaryMassCmp = variable;
    matrixData.boundaryMass = strcat(domain,'_',matrixData.boundaryMassCmp,'_bdM');
    edgename = strcat(domain,'_',matrixData.boundaryMassCmp,'_1_robin');
    [A] = boundaryMass(A,mesh,dof,domain,matrixData.boundaryMass,edgename,variable);
    A.(matrixData.C) = A.(matrixData.C) + laplaceData(24)*A.(matrixData.boundaryMass);  % A
end

% Assemble matrix SDIRK
[a,b,c] = sdirkArray(settings.sdirk);
matrixData.sdirk = strcat(domain,'_M'); % overwrites the mass matrix (never used alone afterwards)
A.(matrixData.sdirk) = A.(matrixData.M) + a(1,1)*dt*A.(matrixData.C);

% Handling Dirichlet boundary conditions
if dof.(bcType)(1)
    dirichletIndices = strcat(domainVariable,'_1_dirichlet');
    dataFlagDirichlet = 21;
    matrixData.dirichlet = strcat(domain,'_dirichlet');
    % coordinated of Dirichlet nodes
    xDirichlet = mesh.(nodes)(1,dof.(dirichletIndices));
    yDirichlet = mesh.(nodes)(2,dof.(dirichletIndices));
    % Extract matrix to impose Dirichlet boundary conditions
    A.(matrixData.dirichlet) = A.(matrixData.C)(dof.(internalIndices),dof.(dirichletIndices));
end

% Select only internal nodes for the SDIRK matrix and factorise
A.(matrixData.sdirk) = A.(matrixData.sdirk)(dof.(internalIndices),dof.(internalIndices));
A.(matrixData.sdirk) = chol(A.(matrixData.sdirk));
% Select only internal rows for the stiffness matrix (used at rhs)
A.(matrixData.C) = A.(matrixData.C)(dof.(internalIndices),dof.(internalIndices));

% Initialise the matrix containing K
s = settings.sdirk;
K = zeros(dof.(numberOfInternalNodes),s);
% Time-advancing loop
for timeStep = 2:N
    % Time instant "tn":
    tn  = t(timeStep-1,1);
    % Numerical solution "un" at time "tn":
    un = u.(domainVariable)(:,timeStep-1);
    un = un(dof.(internalIndices),1);
    % Compute Ki for all the stages
    for i = 1:s
        % time at stage i
        ti = tn + c(i,1)*dt;
        % Assemble rhs
        whichMatrix = [0,0,1,1];
        [A,rhs,matrixData] = laplaceMatrices(mesh,dof,domain,variable,ti,@laplaceData,matrixData,whichMatrix,A);
        rhs.d1_v1 = rhs.d1_v1 - rhs.d1_DtDirichlet;
        % Impose possible Neumann boundary conditions
        if dof.(bcType)(2)
            rhs = laplaceNaturalBC(rhs,mesh,dof,domain,variable,@laplaceData,'neumann',ti);
        end
        % Impose possible Robin boundary conditions
        if dof.(bcType)(3)
            rhs = laplaceNaturalBC(rhs,mesh,dof,domain,variable,@laplaceData,'robin',ti);
        end
        % Restrict the rhs to internal nodes only
        rhs.d1_v1 = rhs.d1_v1(dof.(internalIndices),:);
        % Impose possible Dirichlet boundary conditions
        if dof.(bcType)(1)
            dirichletValuesTi = laplaceData(dataFlagDirichlet,xDirichlet,yDirichlet,ti)';
            rhs.d1_v1 = rhs.d1_v1 - A.(matrixData.dirichlet)*dirichletValuesTi;
        end
        % Sum the contributions from the previous stages
        aK = zeros(dof.(numberOfInternalNodes),1);
        for j = 1:i-1
            aK = aK + a(i,j)*K(:,j);
        end
        % Update "un" at stage "i":
        ui = un + dt*aK;
        % Multiply A*ui
        rhs.d1_v1 = rhs.d1_v1 - A.(matrixData.C)*ui;
        % Compute Ki solving the linear system        
        K(:,i) = (A.(matrixData.sdirk)')\rhs.d1_v1;
        K(:,i) = A.(matrixData.sdirk)\K(:,i);
    end
    % Sum the contributions from the different stages:
    bK = zeros(dof.(numberOfInternalNodes),1);
    for i=1:s
        bK = bK + b(1,i)*K(:,i);
    end
    % Compute the solution "un1" at time step tn+1:
    un1 = un + dt*bK;
    u.(domainVariable)(dof.(internalIndices),timeStep) = un1;
    % Dirichlet values at time t_{n+1}
    dirichletValuesTn1 = laplaceData(dataFlagDirichlet,xDirichlet,yDirichlet,t(timeStep,1))';
    u.(domainVariable)(dof.(dirichletIndices),timeStep) = dirichletValuesTn1;
end

% Plot of the solution
if settings.plot
    tIndex = round(N/2);
    w.(domainVariable) = u.(domainVariable)(:,tIndex);
    [X,Y,U] = preparePlotMeshQ (mesh,domain,w,variable,[]);
    figure;
    sup = surf(X,Y,U);
    colormap(jet); set(sup,'EdgeColor','black','FaceColor','interp');
    grid on; title('Computed solution');
    if settings.plotExact
        time = t0 + (tIndex - 1)*dt;
        U = laplaceData(30,X,Y,time);
        figure;
        sup = surf(X,Y,U);
        colormap(jet); set(sup,'EdgeColor','black','FaceColor','interp');
        grid on; title('Exact pressure');
    end
    % Case of TRIANGULAR elements
    % nodes     = strcat('d',domain,'_p');
    % triangles = strcat('d',domain,'_t');
    % solution1 = strcat(dName,'_v1');
    % figure;
    % sup = pdesurf(mesh.(nodes),mesh.(triangles),u.(solution1));
    % colormap(jet); set(sup,'EdgeColor','black');
    % grid on; title('Computed solution');
end

% Compute error estimates (if required)
if settings.error
    error = laplaceError(u,mesh,domain,variable,@laplaceData);
else
    error = [];
end

end