function [dof] = SetNodesLaplaceStokesSP (dof,domains,strategy)

% We assume that:
%  Domain '1' corresponds to Stokes
%  Domain '2' corresponds to Laplace

switch strategy
    case 'Direct'
        % no modification in the internal nodes
        for i=1:length(domains)
            name1 = strcat('d',domains{i},'_dirichlet_1');
            name2 = strcat('d',domains{i},'_internal_1');
            name3 = strcat('d',domains{i},'_gamma');
            if isfield(dof,name1) % There are Dirichlet b.c.
                name4 = strcat('d',domains{i},'_gamma_1');
                % Nodes on Gamma without Dirichlet nodes
                dof.(name4) = setdiff(dof.(name3),dof.(name1));
            end
            name1 = strcat('d',domains{i},'_ninternal_1');
            dof.(name1) = length(dof.(name2));
            %
            name1 = strcat('d',domains{i},'_dirichlet_2');
            if isfield(dof,name1) % Different conditions on the two components of the velocity (for Stokes)
                name2 = strcat('d',domains{i},'_internal_2');
                name3 = strcat('d',domains{i},'_gamma');
                if isfield(dof,name1) % There are Dirichlet b.c.
                    name4 = strcat('d',domains{i},'_gamma_2');
                    % Nodes on Gamma without Dirichlet nodes
                    dof.(name4) = setdiff(dof.(name3),dof.(name1));
                end
                name1 = strcat('d',domains{i},'_ninternal_2');
                dof.(name1) = length(dof.(name2));
            end
        end
    case 'InterfaceVelocity' %,'SchurPressure','DirichletDirichlet','NeumannNeumann'}
        for i=1:length(domains)
            name1 = strcat('d',domains{i},'_dirichlet_1');
            name2 = strcat('d',domains{i},'_internal_1');
            name3 = strcat('d',domains{i},'_gamma');
            name4 = strcat('d',domains{i},'_internal_noG_1');
            name5 = strcat('d',domains{i},'_gamma_1');
            if isfield(dof,name1) % There are Dirichlet b.c.
                % Nodes on Gamma without Dirichlet nodes
                dof.(name5) = setdiff(dof.(name3),dof.(name1));
            else
                dof.(name5) = dof.(name3);
            end
            % Internal nodes without nodes on Gamma
            dof.(name4) = setdiff(dof.(name2),dof.(name5));
            %
            name1 = strcat('d',domains{i},'_ninternal_1');
            dof.(name1) = length(dof.(name2));
            name1 = strcat('d',domains{i},'_ninternal_noG_1');
            dof.(name1) = length(dof.(name4));
            %
            if strcmp(domains{i},'1') % only for Stokes!
                name1 = strcat('d',domains{i},'_dirichlet_2');
                if isfield(dof,name1) % Different conditions on the two components of the velocity (for Stokes)
                    name2 = strcat('d',domains{i},'_internal_2');
                    name3 = strcat('d',domains{i},'_gamma');
                    name4 = strcat('d',domains{i},'_internal_noG_2');
                    name5 = strcat('d',domains{i},'_gamma_2');
                    if isfield(dof,name1) % There are Dirichlet b.c.
                        % Nodes on Gamma without Dirichlet nodes
                        dof.(name5) = setdiff(dof.(name3),dof.(name1));
                    else
                        dof.(name5) = dof.(name3);
                    end
                    % Internal nodes without nodes on Gamma
                    dof.(name4) = setdiff(dof.(name2),dof.(name5));
                    %
                    name1 = strcat('d',domains{i},'_ninternal_2');
                    dof.(name1) = length(dof.(name2));
                    name1 = strcat('d',domains{i},'_ninternal_noG_2');
                    dof.(name1) = length(dof.(name4));
                end
            end
        end
end

return