function [A] = robinbcMatricesStokes (A,matrixData,mesh,dof,domain,data)

bcType = strcat(domain,'_bctype');
numberOfNodesVelocity = strcat(domain,'_nv1');

% ROBIN boundary conditions
if dof.(bcType)(3)==1
    % There are Robin bc on one of the components
    component = zeros(1,2); m = 0;
    domainVelocity = strcat(domain,'_v1');
    % First component of the velocity
    naturalName = strcat(domainVelocity,'_1_robin');
    for i = 1:size(dof.(naturalName),2)
        if i==1
            m = m+1;
            component(m) = 1;
        end
    end
    % Second component of the velocity
    naturalName = strcat(domainVelocity,'_2_robin');
    for i = 1:size(dof.(naturalName),2)
        if i==1
            m = m+1;
            component(m) = 2;
        end
    end
    if m>0
        component = component(1:m);
    else
        component = [];
    end
    %
    for i = 1:size(component,2)
        matrixData.boundaryMassCmp = strcat('v1_',int2str(component(i)));
        matrixData.boundaryMass = strcat(domain,'_',matrixData.boundaryMassCmp,'_bdM');
        edgename = strcat(domain,'_',matrixData.boundaryMassCmp,'_robin');
        % Modify the stiffness matrix
        if component(i)==1
            indices = 1:dof.(numberOfNodesVelocity);
        elseif component(i)==2
            indices = dof.(numberOfNodesVelocity)+1:2*dof.(numberOfNodesVelocity);
        end
        if ~isfield(A,matrixData.boundaryMass)
            [A] = boundaryMass(A,mesh,dof,domain,matrixData.boundaryMass,edgename,'v1');
        end
        A.(matrixData.C11)(indices,indices) = A.(matrixData.C11)(indices,indices) + ...
            data(24,[],[],[],component(i))*A.(matrixData.boundaryMass);
        A = rmfield(A,matrixData.boundaryMass);
    end
end

end