function [u] = SetDarcyBC (un,dirichlet,dof,domain,varargin)

% SetDarcyBC sets the Dirichlet values to the corresponding degrees of
% freedom. This functions creates a structure U which contains the
% following field:
% U.d?_3 > pressure solution of the Darcy problem

if nargin==5
    u = varargin{1};
end

for i=1:length(domain)
    dname  = strcat('d',domain(i));
    int1   = strcat(dname,'_internal');
    dir1   = strcat(dname,'_dirichlet_v3');
    pnodes = strcat(dname,'_v3');
    name   = strcat(dname,'_v3');
    u.(name) = zeros(dof.(pnodes),1);
    if length(un)==dof.(pnodes)
        u.(name) = un;
    else
        N1 = dof.(int1)(1,1);
        u.(name)(dof.(int1)(1,2:N1)) = un;
    end
    u.(name)(dof.(dir1)) = dirichlet.(name);
end

return