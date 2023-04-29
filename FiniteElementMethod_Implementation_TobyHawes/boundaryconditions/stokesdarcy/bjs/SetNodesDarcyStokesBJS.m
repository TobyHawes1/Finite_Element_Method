function [dof] = SetNodesDarcyStokesBJS (dof,domains,strategy)

% We assume that:
%  Domain '1' corresponds to Stokes
%  Domain '2' corresponds to Darcy

switch strategy
    case 'Direct'
        % no modification in the internal nodes
        for i=1:length(domains)
            if strcmp(domains{i},'1')
                % Stokes
                upflag = '1';
            elseif strcmp(domains{i},'2')
                % Darcy
                upflag = '3';
            end
            name1  = strcat('d',domains{i},'_dirichlet_',upflag);
            name2  = strcat('d',domains{i},'_internal_',upflag); %%%
            name3a = strcat('d',domains{i},'_gamma');
            name3b = strcat('d',domains{3-i},'_gamma');
            name5  = strcat('d',domains{i},'_gamma_',upflag);
            name6  = strcat('d',domains{i},'_ninternal_',upflag);
            %
            dof.(name5)(i,:)   = dof.(name3a);
            dof.(name5)(3-i,:) = dof.(name3b);
            %
            if isfield(dof,name1) % There are Dirichlet b.c.
                % Nodes on Gamma without Dirichlet nodes
                [tmp,indices] = setdiff(dof.(name3a),dof.(name1));
                clear tmp;
                dof.(name5) = dof.(name5)(:,indices);
            end
            % Internal nodes without nodes on Gamma
            dof.(name6) = length(dof.(name2));
            %
            if strcmp(domains{i},'1') % only for Stokes!
                name1 = strcat('d',domains{i},'_dirichlet_2');
                if isfield(dof,name1)
                    % Different conditions on the two components of the 
                    % velocity (for Stokes)
                    name2 = strcat('d',domains{i},'_internal_2');
                    name3a = strcat('d',domains{i},'_gamma');
                    name3b = strcat('d',domains{3-i},'_gamma');
                    name5 = strcat('d',domains{i},'_gamma_2');
                    %
                    dof.(name5)(i,:)   = dof.(name3a);
                    dof.(name5)(3-i,:) = dof.(name3b);
                    %
                    % Nodes on Gamma without Dirichlet nodes
                    [tmp,indices] = setdiff(dof.(name3a),dof.(name1));
                    clear tmp;
                    dof.(name5) = dof.(name5)(:,indices);
                    % Internal nodes without nodes on Gamma
                    name1 = strcat('d',domains{i},'_ninternal_2');
                    dof.(name1) = length(dof.(name2));
                end
            end
        end
    case {'SchurVelocity','DirichletXVF'}
        for i=1:length(domains)
            if strcmp(domains{i},'1')
                % Stokes
                upflag = '1';
            elseif strcmp(domains{i},'2')
                % Darcy
                upflag = '3';
            end
            name1  = strcat('d',domains{i},'_dirichlet_',upflag);
            name2  = strcat('d',domains{i},'_internal_',upflag);
            name3a = strcat('d',domains{i},'_gamma');
            name3b = strcat('d',domains{3-i},'_gamma');
            name5  = strcat('d',domains{i},'_gamma_',upflag);
            name6  = strcat('d',domains{i},'_ninternal_',upflag);
            name4  = strcat('d',domains{i},'_internal_noG_',upflag);
            name7  = strcat('d',domains{i},'_ninternal_noG_',upflag);
            %
            dof.(name5)(i,:)   = dof.(name3a);
            dof.(name5)(3-i,:) = dof.(name3b);
            %
            if isfield(dof,name1) % There are Dirichlet b.c.
                % Nodes on Gamma without Dirichlet nodes
                [tmp,indices] = setdiff(dof.(name3a),dof.(name1));
                clear tmp;
                dof.(name5) = dof.(name5)(:,indices);
            end
            % Internal nodes without nodes on Gamma
            dof.(name4) = setdiff(dof.(name2),dof.(name5)(i,:));
            dof.(name6) = length(dof.(name2));
            dof.(name7) = length(dof.(name4));
            %
            if strcmp(domains{i},'1') % only for Stokes!
                name1 = strcat('d',domains{i},'_dirichlet_2');
                if isfield(dof,name1)
                    % Different conditions on the two components of the 
                    % velocity (for Stokes)
                    name2 = strcat('d',domains{i},'_internal_2');
                    name3a = strcat('d',domains{i},'_gamma');
                    name3b = strcat('d',domains{3-i},'_gamma');
                    name4 = strcat('d',domains{i},'_internal_noG_2');
                    name5 = strcat('d',domains{i},'_gamma_2');
                    %
                    dof.(name5)(i,:)   = dof.(name3a);
                    dof.(name5)(3-i,:) = dof.(name3b);
                    %
                    % Nodes on Gamma without Dirichlet nodes
                    [tmp,indices] = setdiff(dof.(name3a),dof.(name1));
                    clear tmp;
                    dof.(name5) = dof.(name5)(:,indices);
                    % Internal nodes without nodes on Gamma
                    dof.(name4) = setdiff(dof.(name2),dof.(name5)(i,:));
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