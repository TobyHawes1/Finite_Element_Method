function [rhs] = darcyScalarStokesBjsCouplingBC (A,rhs,dof,dirichlet,domains)

stokesDomain = strcat('d',domains{1});
stokesRhs    = strcat(stokesDomain,'_v1_v2');
stokesValues = strcat(stokesDomain,'_v2');
darcyDomain  = strcat('d',domains{2});
darcyRhs     = strcat(darcyDomain,'_v3');
darcyValues  = strcat(darcyDomain,'_v3');

dirichletNodesDarcy  = strcat(darcyDomain,'_dirichlet_v3');
dirichletNodesStokes = strcat(stokesDomain,'_dirichlet_v2');
if ~isfield(dof,dirichletNodesStokes)
    dirichletNodesStokes = strcat(stokesDomain,'_dirichlet_v1');
end

unodesS = strcat(stokesDomain,'_v1');

% Correction to the Stokes problem due to the coupling
internalNodesStokes = strcat(stokesDomain,'_internal');
N1 = dof.(internalNodesStokes)(1,1);
N2 = dof.(internalNodesStokes)(1,2);
rhs.(stokesRhs)(dof.(internalNodesStokes)(1,N1+1:N2),1) = ...
    rhs.(stokesRhs)(dof.(internalNodesStokes)(1,N1+1:N2),1) + ...
    A.MGammaSD(dof.(internalNodesStokes)(1,N1+1:N2)-dof.(unodesS),...
                       dof.(dirichletNodesDarcy))*dirichlet.(darcyValues);

% Correction to the Darcy problem due to the coupling
internalNodesDarcy = strcat(darcyDomain,'_internal');
N1 = dof.(internalNodesDarcy)(1,1);
rhs.(darcyRhs)(dof.(internalNodesDarcy)(1,2:N1),1) = ...
    rhs.(darcyRhs)(dof.(internalNodesDarcy)(1,2:N1),1) - ...
    ((A.MGammaSD(dof.(dirichletNodesStokes),dof.(internalNodesDarcy)(1,2:N1)))')*...
    dirichlet.(stokesValues);

end