function [A,rhs,dirichlet] = ImposeDarcyBC (A,matdata,rhs,mesh,dof,domain,...
    data,dirichlet,varargin)

dname  = strcat('d',domain);
bctype = strcat(dname,'_bctype');

% NEUMANN boundary conditions
if dof.(bctype)(2)
    rhs = DarcyNaturalBC(rhs,mesh,dof,domain,data,'neumann',varargin{:});
end

% ROBIN boundary conditions
if dof.(bctype)(3)
    rhs = DarcyNaturalBC(rhs,mesh,dof,domain,data,'robin',varargin{:});
    % Modify the stiffness matrix
    matdata.BoundaryMassCmp = '1';
    matdata.BoundaryMass = strcat(dname,'_bdM_',matdata.BoundaryMassCmp);
    edgename = strcat(dname,'_robin_',matdata.BoundaryMassCmp);
    if ~isfield(A,matdata.BoundaryMass)
        [A] = BoundaryMass(A,mesh,dof,domain,matdata.BoundaryMass,edgename);
    end
    A.(matdata.C) = A.(matdata.C) + data(24)*A.(matdata.BoundaryMass);
end

% DIRICHLET boundary conditions
if dof.(bctype)(1)
    [rhs,dirichlet] = DarcyDirichletBC(A,matdata,rhs,mesh,dof,domain,data,dirichlet);
end

if (dof.(bctype)(1)==0 & dof.(bctype)(3)==0)
    fprintf('   [Darcy problem] Warning: solution defined up to an additive constant!\n');
end

return