function [rhs,dirichlet] = stokesDirichletBC (A,matrixData,rhs,mesh,dof,domain,...
    data,settings,dirichlet,time)

% StokesDirichletBC imposes Dirichlet boundary conditions for the Stokes
% problem.

domainVelocity = strcat(domain,'_v1');
nodesVelocity  = strcat(domainVelocity,'_p');
dirichlet1 = strcat(domainVelocity,'_1_dirichlet');
dirichlet2 = strcat(domainVelocity,'_2_dirichlet');
numberOfNodesVelocity = strcat(domain,'_nv1');
diri1 = strcat(domain,'_v1_1');
diri2 = strcat(domain,'_v1_2');
internalNodes  = strcat(domainVelocity,'_internal');
domainPressure = strcat(domain,'_v2');

% internalName   = strcat(dName,'_internal');
% dirichletName2 = strcat(dName,'_dirichlet_v2');
% if ~isfield(dof,dirichletName2)
%     dirichletName2 = dirichletName1;
% end
% uNodes = strcat(dName,'_v1');
% nodes  = strcat(dName,'_p');

% Data "21" corresponds to Dirichlet bc
dataFlag = 21;
x = mesh.(nodesVelocity)(1,dof.(dirichlet1));
y = mesh.(nodesVelocity)(2,dof.(dirichlet1));
dirichlet.(diri1) = feval(data,dataFlag,x,y,time,1)';
x = mesh.(nodesVelocity)(1,dof.(dirichlet2));
y = mesh.(nodesVelocity)(2,dof.(dirichlet2));
dirichlet.(diri2) = feval(data,dataFlag,x,y,time,2)';

dirichletValues = [dirichlet.(diri1); dirichlet.(diri2)];
dirichletNodes  = [dof.(dirichlet1) dof.(dirichlet2) + dof.(numberOfNodesVelocity)]';

% Velocity (momentum equation)
rhs.(domainVelocity)(dof.(internalNodes)) = ...
    rhs.(domainVelocity)(dof.(internalNodes)) - ...
    A.(matrixData.C11)(dof.(internalNodes),dirichletNodes)*dirichletValues;
% Velocity (continuity equation)
if settings.stokesStab == 4
    % Douglas-Wang stabilization (non-symmetric linear system)
    rhs.(domainPressure) = rhs.(domainPressure) - ...
        A.(matrixData.D)(:,dirichletNodes)*dirichletValues;
else
    rhs.(domainPressure) = rhs.(domainPressure) - ...
        (A.(matrixData.G)(dirichletNodes,:)')*dirichletValues;
end

end