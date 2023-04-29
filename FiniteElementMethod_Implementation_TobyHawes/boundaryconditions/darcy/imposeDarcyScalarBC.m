function [A,rhs,dof,dirichlet] = imposeDarcyScalarBC (A,matrixData,rhs,mesh,dof,domain,...
    data,dirichlet,varargin)

dname  = strcat('d',domain);
bctype = strcat(dname,'_bctype');

% NEUMANN boundary conditions
if dof.(bctype)(2)==1
    rhs = darcyScalarNaturalBC(rhs,mesh,dof,domain,data,'neumann',varargin{:});
end

% ROBIN boundary conditions
if dof.(bctype)(3)==1
    rhs = darcyScalarNaturalBC(rhs,mesh,dof,domain,data,'robin',varargin{:});
    % Modify the stiffness matrix
    matrixData.boundaryMassCmp = '3';
    matrixData.boundaryMass = strcat(dname,'_bdM_v',matrixData.boundaryMassCmp);
    edgename = strcat(dname,'_robin_v',matrixData.boundaryMassCmp);
    if ~isfield(A,matrixData.boundaryMass)
        [A] = boundaryMass(A,mesh,dof,domain,matrixData.boundaryMass,edgename,'v3');
    end
    A.(matrixData.C) = A.(matrixData.C) + data(24)*A.(matrixData.BoundaryMass);
    A = rmfield(A,matrixData.boundaryMass);
end

% DIRICHLET boundary conditions
if dof.(bctype)(1)==1
    [rhs,dirichlet] = darcyScalarDirichletBC(A,matrixData,rhs,mesh,dof,domain,data,dirichlet);
end

if (dof.(bctype)(1)==0 & dof.(bctype)(3)==0)
    fprintf('   [Darcy problem] Warning: solution defined up to an additive constant!\n');
end

% Clean unused variables
nameN1 = strcat(dname,'_neumann_v3');
nameR1 = strcat(dname,'_robin_v3');
name = {nameN1,nameR1};
for i=1:length(name)
    if isfield(dof,name{i})
        dof = rmfield(dof,name{i});
    end
end

end