function [A,rhs,dof,dirichlet] = imposeStokesBC (A,matrixData,rhs,...
    mesh,dof,domain,data,dirichlet,settings,time,varargin)

% ImposeStokesBC imposes boundary conditions for the Stokes problem.
% Boundary conditions may be of Dirichlet, Neumann and Robin type.
% Supported finite elements: P1-P1, MINI, Taylor-Hood, and P2-P2.

bcType = strcat(domain,'_bctype');
numberOfNodesVelocity = strcat(domain,'_nv1');

% NEUMANN boundary conditions
if dof.(bcType)(2)==1
    % There are Neumann bc on one of the components
    if ~isempty(varargin)
        arguments = varargin{1};
        rhs = stokesNaturalBC(rhs,mesh,dof,domain,data,'neumann',time,arguments);
    else
        rhs = stokesNaturalBC(rhs,mesh,dof,domain,data,'neumann',time);
    end
end

% ROBIN boundary conditions
if dof.(bcType)(3)==1
    % There are Robin bc on one of the components
    [rhs,component] = stokesNaturalBC(rhs,mesh,dof,domain,data,'robin',time);
    for i = 1:size(component,2)
        matrixData.boundaryMassCmp = strcat('v1_',int2str(component(i)));
        matrixData.boundaryMass = strcat(domain,'_',matrixData.boundaryMassCmp,'_bdM');
        edgename = strcat(domain,'_',matrixData.boundaryMassCmp,'_robin');
        % Modify the stiffness matrix
        if component(i)==1
            indices = 1:dof.(numberOfNodesVelocity);
        elseif component(i)==2
            indices = dof.(numberOfNodesVelocity)+1:2*dof.(numberOfNodesVelocity);
            %if ~isfield(dof,edgename)
            %    edgename = strcat(dName,'_robin_v',int2str(component(1)));
            %end
        end
        if settings.robin==0
            if ~isfield(A,matrixData.boundaryMass)
                [A] = boundaryMass(A,mesh,dof,domain,matrixData.boundaryMass,edgename,'v1');
            end
            A.(matrixData.C11)(indices,indices) = A.(matrixData.C11)(indices,indices) + ...
                data(24,[],[],[],component(i))*A.(matrixData.boundaryMass);
        elseif settings.robin==1
            % THIS NEEDS DEBUGGING!!!
            [A] = boundaryMassNI(A,mesh,dof,domain,...
                matrixData.boundaryMass,edgename,data,component(i));
            A.(matrixData.C11)(indices,indices) = A.(matrixData.C11)(indices,indices) + ...
                A.(matrixData.boundaryMass);
        end
        A = rmfield(A,matrixData.boundaryMass);
    end
end

% DIRICHLET boundary conditions
if dof.(bcType)(1)==1
    % There are Dirichlet bc
    [rhs,dirichlet] = stokesDirichletBC(A,matrixData,rhs,mesh,dof,domain,data,...
        settings,dirichlet,time);
else
    dirichlet = {};
end

% % Clear unused variables
% domainVelocity = strcat(domain,'_v1');
% nameN1 = strcat(domainVelocity,'_1_neumann');
% nameN2 = strcat(domainVelocity,'_2_neumann');
% nameR1 = strcat(domainVelocity,'_1_robin');
% nameR2 = strcat(domainVelocity,'_2_robin');
% name = {nameN1,nameN2,nameR1};%,nameR2};
% for i=1:length(name)
%     if isfield(dof,name{i})
%         dof = rmfield(dof,name{i});
%     end
% end

end