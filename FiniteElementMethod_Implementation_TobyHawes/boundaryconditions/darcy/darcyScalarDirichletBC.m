function [rhs,dirichlet] = darcyScalarDirichletBC (A,matrixData,rhs,mesh,dof,...
    domain,data,dirichlet,varargin)

dName         = strcat('d',domain);
name3         = strcat(dName,'_v3');
dirichletName = strcat(dName,'_dirichlet_v3');
internalD     = strcat(dName,'_internal');
nodes         = strcat(dName,'_p');
diri1         = strcat(dName,'_v3');

% Data "22" corresponds to Dirichlet bc
dataflag = 22;
x = mesh.(nodes)(1,dof.(dirichletName));
y = mesh.(nodes)(2,dof.(dirichletName));
dirichlet.(diri1) = feval(data,dataflag,x,y,varargin{:})';

N = dof.(internalD)(1,1);
rhs.(name3)(dof.(internalD)(1,2:N)) = rhs.(name3)(dof.(internalD)(1,2:N)) - ...
    A.(matrixData.C)(dof.(internalD)(1,2:N),dof.(dirichletName))*dirichlet.(diri1);

end