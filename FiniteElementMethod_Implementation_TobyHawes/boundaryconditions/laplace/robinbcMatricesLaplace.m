function [A] = robinbcMatricesLaplace (A,matrixData,mesh,dof,domain,variable,data,varargin)

domainVariable = strcat(domain,'_',variable);
bcType = strcat(domainVariable,'_bctype');

% ROBIN boundary conditions
if dof.(bcType)(3)
    %rhs = laplaceNaturalBC(rhs,mesh,dof,domain,variable,data,'robin',time,varargin{:});
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

end