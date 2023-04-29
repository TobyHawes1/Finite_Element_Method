function [u] = setStokesBC (un,dirichlet,dof,domain,varargin)

% SETSTOKESBC sets the Dirichlet values to the corresponding degrees of
% freedom. This functions creates a structure U which contains the
% following fields:
% U.d?_1 > first component of the velocity
% U.d?_2 > second component of the velocity
% U.d?_3 > pressure

if ~isempty(varargin)
    u = varargin{1};
end

domainVelocity = strcat(domain,'_v1');
dirichlet1  = strcat(domainVelocity,'_1_dirichlet');
dirichlet2  = strcat(domainVelocity,'_2_dirichlet');
dirichletValues1 = strcat(domainVelocity,'_1');
dirichletValues2 = strcat(domainVelocity,'_2');
numberOfNodesVelocity = strcat(domain,'_nv1');
numberOfNodesPressure = strcat(domain,'_nv2');
velocity1   = strcat(domain,'_v1_1');
velocity2   = strcat(domain,'_v1_2');
pressure    = strcat(domain,'_v2');

if length(un) == 2*dof.(numberOfNodesVelocity)+dof.(numberOfNodesPressure)
    u.(velocity1) = un(1:dof.(numberOfNodesVelocity));
    u.(velocity2) = un(dof.(numberOfNodesVelocity)+1:2*dof.(numberOfNodesVelocity));
    u.(pressure)  = un(2*dof.(numberOfNodesVelocity)+1:end);
else
    internalNodes  = strcat(domainVelocity,'_internal');
    numberOfInternalNodes = strcat(domainVelocity,'_ninternal');
    % Velocity
    u.(velocity1) = zeros(2*dof.(numberOfNodesVelocity),1);
    u.(velocity1)(dof.(internalNodes)) = un(1:dof.(numberOfInternalNodes));
    u.(velocity2) = u.(velocity1)(dof.(numberOfNodesVelocity)+1:end);
    u.(velocity1) = u.(velocity1)(1:dof.(numberOfNodesVelocity));
    % Pressure
    u.(pressure) = un(dof.(numberOfInternalNodes)+1:dof.(numberOfInternalNodes)+dof.(numberOfNodesPressure));
end
clear un;
% Set Dirichlet values
if ( isfield(dof,dirichlet1) && ~isempty(dof.(dirichlet1)) )
    u.(velocity1)(dof.(dirichlet1)) = dirichlet.(dirichletValues1);
end
if ( isfield(dof,dirichlet2) && ~isempty(dof.(dirichlet2)) )
    u.(velocity2)(dof.(dirichlet2)) = dirichlet.(dirichletValues2);
end

end