function [u] = setLaplaceBC (un,dirichlet,dof,domain,variable,varargin)

% SETLAPLACEBC sets the Dirichlet values to the corresponding degrees of
% freedom. This functions creates a structure U which contains the
% following field:
% U.d?_1 > solution of the Laplace problem

if ~isempty(varargin)
    u = varargin{1};
end

domainVariable = strcat(domain,'_',variable);

int1   = strcat(domainVariable,'_1_internal');
dir1   = strcat(domainVariable,'_1_dirichlet');

numberOfNodes = strcat(domain,'_n',variable);

u.(domainVariable) = zeros(dof.(numberOfNodes),1);
if length(un)==dof.(domainVariable)
    u.(domainVariable) = un;
else
    u.(domainVariable)(dof.(int1)) = un;
end
u.(domainVariable)(dof.(dir1)) = dirichlet.(domainVariable);

end