function [rhs] = stokesdarcybjsCouplingBC (A,rhs,dof,dirichlet,domainStokes,domainDarcy)

rhsStokes    = strcat(domainStokes,'_v1');
stokesValues = strcat(domainStokes,'_v1_2');
rhsDarcy     = strcat(domainDarcy,'_v1');
darcyValues  = strcat(domainDarcy,'_v1');

dirichletNodesStokes2 = strcat(domainStokes,'_v1_2_dirichlet');
dirichletNodesDarcy   = strcat(domainDarcy,'_v1_1_dirichlet');

% Correction to the Stokes problem due to the coupling
numberOfInternalNodesStokesVelocity1 = strcat(domainStokes,'_v1_1_ninternal');
internalNodesStokes  = strcat(domainStokes,'_v1_internal');
internalNodesStokes2 = strcat(domainStokes,'_v1_2_internal');
indicesStokes2 = dof.(internalNodesStokes)(1,dof.(numberOfInternalNodesStokesVelocity1)+1:end);
rhs.(rhsStokes)(indicesStokes2,1) = rhs.(rhsStokes)(indicesStokes2,1) + ...
    A.MGammaSD(dof.(internalNodesStokes2),dof.(dirichletNodesDarcy))*dirichlet.(darcyValues);

% Correction to the Darcy problem due to the coupling
internalNodesDarcy = strcat(domainDarcy,'_v1_1_internal');
rhs.(rhsDarcy)(dof.(internalNodesDarcy),1) = ...
    rhs.(rhsDarcy)(dof.(internalNodesDarcy),1) - ...
    ((A.MGammaSD(dof.(dirichletNodesStokes2),dof.(internalNodesDarcy)))')*...
    dirichlet.(stokesValues);

end