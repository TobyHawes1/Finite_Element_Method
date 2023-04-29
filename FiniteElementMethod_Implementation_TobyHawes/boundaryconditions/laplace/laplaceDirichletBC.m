function [rhs,dirichlet] = laplaceDirichletBC (A,matrixData,rhs,mesh,dof,domain,variable,time,data,dirichlet,varargin)

domainVariable = strcat(domain,'_',variable);
dirichletName = strcat(domainVariable,'_1_dirichlet');
internalName  = strcat(domainVariable,'_1_internal');
nodes = strcat(domainVariable,'_p');

x = mesh.(nodes)(1,dof.(dirichletName));
y = mesh.(nodes)(2,dof.(dirichletName));
% dataFlag 21 corresponds to Dirichlet bc
dataFlagDirichlet= 21;
dirichlet.(domainVariable) = feval(data,dataFlagDirichlet,x,y,time,varargin{:})';

rhs.(domainVariable)(dof.(internalName)) = rhs.(domainVariable)(dof.(internalName)) - ...
    A.(matrixData.C)(dof.(internalName),dof.(dirichletName))*dirichlet.(domainVariable);

end