function [rhs,dirichlet] = DarcyDirichletBC (A,matdata,rhs,mesh,dof,domain,data,dirichlet,varargin)

dname = strcat('d',domain);
dirichlet_name1 = strcat(dname,'_dirichlet_3');
internal_name1  = strcat(dname,'_internal_3');
nodes = strcat(dname,'_p');
diri1 = strcat(dname,'_3');

x = mesh.(nodes)(1,dof.(dirichlet_name1));
y = mesh.(nodes)(2,dof.(dirichlet_name1));
% Data "21" corresponds to Dirichlet bc
dirichlet.(diri1) = feval(data,21,x,y,varargin{:})';

rhs.(dname)(dof.(internal_name1)) = rhs.(dname)(dof.(internal_name1)) - ...
    A.(matdata.C)(dof.(internal_name1),dof.(dirichlet_name1))*dirichlet.(diri1);

return