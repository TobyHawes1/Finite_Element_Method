function [u] = SetDarcyScalarBC (un,dirichlet,dof,domain,varargin)

% SetDarcyScalarBC sets the Dirichlet values to the corresponding degrees
% of freedom. This functions creates a structure U which contains the
% following field:
% U.d?_3 > pressure solution of the Darcy problem

if nargin==5
    u = varargin{1};
end

for i=1:length(domain)
    dname      = strcat('d',domain(i));
    int        = strcat(dname,'_internal');
    dirichlet1 = strcat(dname,'_dirichlet_v3');
    pNodes     = strcat(dname,'_v3');
    name       = strcat(dname,'_v3');
    N          = dof.(int)(1,1);
    u.(name)   = zeros(dof.(pNodes),1);
    if length(un)==dof.(pNodes)
        u.(name) = un;
    else
        u.(name)(dof.(int)(1,2:N)) = un;
    end
    u.(name)(dof.(dirichlet1)) = dirichlet.(name);
end

end