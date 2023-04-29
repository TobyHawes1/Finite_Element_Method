function [A,rhs,dirichletValues] = imposeLaplaceBC (A,matrixData,rhs,mesh,dof,domain,variable,...
    time,data,dirichletValues,varargin)

domainVariable = strcat(domain,'_',variable);
bcType = strcat(domainVariable,'_bctype');

% NEUMANN boundary conditions
if dof.(bcType)(2)
    rhs = laplaceNaturalBC(rhs,mesh,dof,domain,variable,data,'neumann',time,varargin{:});
end

% ROBIN boundary conditions
if dof.(bcType)(3)
    rhs = laplaceNaturalBC(rhs,mesh,dof,domain,variable,data,'robin',time,varargin{:});
    % Modify the stiffness matrix
    matrixData.boundaryMassCmp = variable;
    matrixData.boundaryMass = strcat(domain,'_',matrixData.boundaryMassCmp,'_bdM');
    edgename = strcat(domain,'_',matrixData.boundaryMassCmp,'_1_robin');
    if ~isfield(A,matrixData.boundaryMass)
        [A] = boundaryMass(A,mesh,dof,domain,matrixData.boundaryMass,edgename,variable);
    end
    A.(matrixData.C) = A.(matrixData.C) + data(24)*A.(matrixData.boundaryMass);
    A = rmfield(A,matrixData.boundaryMass);
end

% DIRICHLET boundary conditions
if dof.(bcType)(1)
    [rhs,dirichletValues] = laplaceDirichletBC(A,matrixData,rhs,mesh,dof,domain,variable,time,data,dirichletValues);
end

% if (dof.(bcType)(1)==0 && dof.(bcType)(3)==0)
%     fprintf('   [Laplace problem] Warning: solution defined up to an additive constant!\n');
% end

% Clean unused variables
% nameN1 = strcat(domainVariable,'_1_neumann');
% nameR1 = strcat(domainVariable,'_1_robin');
% name = {nameN1,nameR1};
% for i=1:length(name)
%     if isfield(dof,name{i})
%         dof = rmfield(dof,name{i});
%     end
% end

end