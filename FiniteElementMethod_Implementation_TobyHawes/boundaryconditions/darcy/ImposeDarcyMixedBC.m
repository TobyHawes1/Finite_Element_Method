function [A,rhs,dof,dirichlet] = ImposeDarcyMixedBC (A,matdata,rhs,...
    mesh,dof,domain,data,dirichlet,settings,varargin)

% ImposeDarcyMixedBC imposes boundary conditions for the Stokes problem.
% Boundary conditions may be of Dirichlet or Neumann type.
% Supported finite elements: P2-P1 (Taylor-Hood) and P1b-P1 (MINI).

dname  = strcat('d',domain);
bctype = strcat(dname,'_bctype');

% NEUMANN boundary conditions
if dof.(bctype)(2)
    % There are Neumann bc on one of the components
    rhs = DarcyMixedNaturalBC(rhs,mesh,dof,domain,data,'neumann',settings,...
        varargin{:});
end

% DIRICHLET boundary conditions
if dof.(bctype)(1)
    % There are Dirichlet bc
    [rhs,dirichlet,dof] = DarcyMixedDirichletBC(A,matdata,rhs,...
        mesh,dof,domain,data,settings,dirichlet,varargin{:});
end

return