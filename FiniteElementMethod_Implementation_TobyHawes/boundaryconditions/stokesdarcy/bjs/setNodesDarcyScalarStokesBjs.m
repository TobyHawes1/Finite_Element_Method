function [dof,settings] = setNodesDarcyScalarStokesBjs (dof,domains,settings)

% We assume that:
%  Domain '1' corresponds to Stokes
%  Domain '2' corresponds to Darcy

switch settings.method
    case {'SchurVelocity','SchurPressure'}
        % ----------------------------------------------------------------
        % STOKES PROBLEM
        % ----------------------------------------------------------------
        dName = strcat('d',domains{1});
        % Velocity (second component) interface nodes without Dirichlet
        % nodes
        dirichletNodes = strcat(dName,'_dirichlet_v2');
        if ~isfield(dof,dirichletNodes)
            dirichletNodes = strcat(dName,'_dirichlet_v1');
        end
        gammaNodes = strcat(dName,'_gamma_v1');
        dof.(gammaNodes) = setdiff(dof.(gammaNodes),dof.(dirichletNodes));
        % Internal nodes without interface nodes for the second component
        % of the velocity
        uNodes = strcat(dName,'_v1');
        internalNodes = strcat(dName,'_internal');
        N1 = dof.(internalNodes)(1,1);
        N2 = dof.(internalNodes)(1,2);
        L1 = length(dof.(gammaNodes));
        nint2 = N2-N1;
        nint2 = nint2-L1;
        dof.(internalNodes)(1,N1+1:nint2+N1) = ...
            setdiff(dof.(internalNodes)(1,N1+1:N2),...
            dof.(gammaNodes)+dof.(uNodes));
        dof.(internalNodes)(1,2) = nint2+N1;
        N2 = dof.(internalNodes)(1,2);
        dof.(internalNodes) = dof.(internalNodes)(1,1:N2);
        % Remove unused dofs
        %GammaNodes = strcat(dname,'_gamma_v3');
        %dof = rmfield(dof,GammaNodes);
        %
        % ----------------------------------------------------------------
        % DARCY PROBLEM
        % ----------------------------------------------------------------
        if strcmp(settings.method,'SchurVelocity')
            dName = strcat('d',domains{2});
            % Remove unused dofs
            %GammaNodes12 = strcat(dname,'_gamma_v1');
            gammaNodes3  = strcat(dName,'_gamma_v3');
            %dof = rmfield(dof,{GammaNodes12,GammaNodes3});
            dof = rmfield(dof,gammaNodes3);
        else
            dName = strcat('d',domains{2});
            % Pressure interface nodes without Dirichlet nodes
            dirichletNodes = strcat(dName,'_dirichlet_v3');
            gammaNodes = strcat(dName,'_gamma_v3');
            dof.(gammaNodes) = setdiff(dof.(gammaNodes),dof.(dirichletNodes));
            % Internal nodes without interface nodes for the pressure
            internalNodes = strcat(dName,'_internal');
            N1 = dof.(internalNodes)(1,1);
            L1 = length(dof.(gammaNodes));
            nInt = (N1-1)-L1;
            dof.(internalNodes)(1,2:nInt+1) = ...
                setdiff(dof.(internalNodes)(1,2:N1),dof.(gammaNodes));
            dof.(internalNodes)(1,1) = nInt+1;
            N1 = dof.(internalNodes)(1,1);
            dof.(internalNodes) = dof.(internalNodes)(1,1:N1);
        end
        %
%     case {'DirichletXVF'}
%         if strcmp(settings.preconditioner,'OS')
%             % Check if Dirichlet nodes coincide with interface nodes
%             % Stokes domain
%             dname = strcat('d',domains{1});
%             DirichletNodes = strcat(dname,'_dirichlet_2');
%             if ~isfield(dof,DirichletNodes)
%                 DirichletNodes = strcat(dname,'_dirichlet_1');
%             end
%             GammaNodes = strcat(dname,'_gamma_1_2');
%             Intersection = intersect(dof.(DirichletNodes),dof.(GammaNodes));
%             if ~isempty(Intersection)
%                 IntersectionName = strcat(dname,'_corners');
%                 dof.(IntersectionName) = Intersection;
%                 settings.BCCorrection(1) = 1;
%             else
%                 settings.BCCorrection(1) = 0;
%             end
%             % Darcy domain
%             dname = strcat('d',domains{2});
%             DirichletNodes = strcat(dname,'_dirichlet_3');
%             GammaNodes     = strcat(dname,'_gamma_3');
%             Intersection = intersect(dof.(DirichletNodes),dof.(GammaNodes));
%             if ~isempty(Intersection)
%                 IntersectionName = strcat(dname,'_corners');
%                 dof.(IntersectionName) = Intersection;
%                 settings.BCCorrection(2) = 1;
%             else
%                 settings.BCCorrection(2) = 0;
%             end
%         end
%         % STOKES PROBLEM
%         % ----------------------------------------------------------------
%         dname = strcat('d',domains{1});
%         % Velocity (second component) interface nodes without Dirichlet
%         % nodes
%         DirichletNodes = strcat(dname,'_dirichlet_2');
%         if ~isfield(dof,DirichletNodes)
%             DirichletNodes = strcat(dname,'_dirichlet_1');
%         end
%         GammaNodes = strcat(dname,'_gamma_1_2');
%         dof.(GammaNodes) = setdiff(dof.(GammaNodes),dof.(DirichletNodes));
%         % Internal nodes without interface nodes for the second component
%         % of the velocity
%         unodes = strcat(dname,'_u');
%         InternalNodes = strcat(dname,'_internal');
%         N1 = dof.(InternalNodes)(1,1);
%         N2 = dof.(InternalNodes)(1,2);
%         L1 = length(dof.(GammaNodes));
%         nint2 = N2-N1;
%         nint2 = nint2-L1;
%         dof.(InternalNodes)(1,N1+1:nint2+N1) = ...
%             setdiff(dof.(InternalNodes)(1,N1+1:N2),...
%             dof.(GammaNodes)+dof.(unodes));
%         dof.(InternalNodes)(1,2) = nint2+N1;
%         N2 = dof.(InternalNodes)(1,2);
%         dof.(InternalNodes) = dof.(InternalNodes)(1,1:N2);
%         % Remove unused dofs
%         GammaNodes = strcat(dname,'_gamma_3');
%         dof = rmfield(dof,GammaNodes);
%         %
%         % DARCY PROBLEM
%         % ----------------------------------------------------------------
%         dname = strcat('d',domains{2});
%         % Pressure interface nodes without Dirichlet nodes
%         DirichletNodes = strcat(dname,'_dirichlet_3');
%         GammaNodes = strcat(dname,'_gamma_3');
%         dof.(GammaNodes) = setdiff(dof.(GammaNodes),dof.(DirichletNodes));
%         % Internal nodes without interface nodes for the pressure
%         InternalNodes = strcat(dname,'_internal');
%         N1 = dof.(InternalNodes)(1,1);
%         L1 = length(dof.(GammaNodes));
%         nInt = (N1-1)-L1;
%         dof.(InternalNodes)(1,2:nInt+1) = ...
%             setdiff(dof.(InternalNodes)(1,2:N1),dof.(GammaNodes));
%         dof.(InternalNodes)(1,1) = nInt+1;
%         N1 = dof.(InternalNodes)(1,1);
%         dof.(InternalNodes) = dof.(InternalNodes)(1,1:N1);
%         % Remove unused dofs
%         GammaNodes = strcat(dname,'_gamma_1_2');
%         dof = rmfield(dof,GammaNodes);
%         %
%         % Obtain correspondence between interface nodes
%         % ----------------------------------------------------------------
%         d1name  = strcat('d',domains{1});
%         p1name = strcat(d1name,'_p');
%         d2name  = strcat('d',domains{2});
%         p2name = strcat(d2name,'_p');
%         % Stokes interface nodes corresponding to the Darcy interface ones
%         % (used for hat{SigmaS})
%         DarcyGammaNodes = strcat(d2name,'_gamma_3');
%         dof.GammaNodesDS = 0*dof.(DarcyGammaNodes);
%         for i=1:length(dof.GammaNodesDS)
%             xyDarcy = mesh.(p2name)(:,dof.(DarcyGammaNodes)(1,i));
%             flag = 0;
%             j = 0;
%             while flag==0
%                 j = j+1;
%                 xyCheck = mesh.(p1name)(:,j);
%                 xyCheck = norm(xyCheck-xyDarcy);
%                 if xyCheck<1.e-7
%                     flag = 1;
%                 end
%             end
%             dof.GammaNodesDS(i) = j;
%         end
%         % Darcy interface nodes corresponding to the Stokes interface ones
%         % (used for hat{SigmaP})
%         StokesGammaNodes = strcat(d1name,'_gamma_1_2');
%         dof.GammaNodesSD = 0*dof.(StokesGammaNodes);
%         for i=1:length(dof.GammaNodesSD)
%             xyStokes = mesh.(p1name)(:,dof.(StokesGammaNodes)(1,i));
%             flag = 0;
%             j = 0;
%             while flag==0
%                 j = j+1;
%                 xyCheck = mesh.(p2name)(:,j);
%                 xyCheck = norm(xyCheck-xyStokes);
%                 if xyCheck<1.e-7
%                     flag = 1;
%                 end
%             end
%             dof.GammaNodesSD(i) = j;
%         end
end

end