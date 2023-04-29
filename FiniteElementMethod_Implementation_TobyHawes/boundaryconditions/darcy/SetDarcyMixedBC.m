function [u] = SetDarcyMixedBC (un,dirichlet,dof,domain,settings,varargin)

% SetDarcyMixedBC sets the Dirichlet values to the corresponding degrees of
% freedom. This functions creates a structure U which contains the
% following fields:
% U.d?_1 > first component of the velocity
% U.d?_2 > second component of the velocity
% U.d?_3 > pressure

if nargin==6
    u = varargin{1};
end

if settings.darcymixed<4
    % Essential boundary conditions are on the (normal) velocity
    for i=1:length(domain)
        dname         = strcat('d',domain(i));
        internal_name = strcat(dname,'_internal');
        %nint1   = strcat(dname,'_ninternal_1');
        dir1    = strcat(dname,'_dirichlet_1');
        dir2    = strcat(dname,'_dirichlet_2');
        ddir1   = strcat(dname,'_1');
        ddir2   = strcat(dname,'_2');
        %pnodes  = strcat(dname,'_all_p');
        nunodes = strcat(dname,'_u');
        npnodes = strcat(dname,'_p');
        vel1    = strcat(dname,'_1');
        vel2    = strcat(dname,'_2');
        press   = strcat(dname,'_3');
        % Allocate memory space
        u.(vel1)  = zeros(2*dof.(nunodes),1);
        u.(vel2)  = zeros(dof.(nunodes),1);
        u.(press) = zeros(dof.(npnodes),1);
        %
        % Velocity
        u.(vel1)(dof.(internal_name)(3:end)) = un(1:dof.(internal_name)(2)-2);
        u.(vel2) = u.(vel1)(dof.(nunodes)+1:end);
        u.(vel1) = u.(vel1)(1:dof.(nunodes));
        % Pressure
        u.(press) = un(dof.(internal_name)(2)-2+1:dof.(internal_name)(2)-2+dof.(npnodes));
        clear un;
        % Set Dirichlet values
        if ( isfield(dof,dir1) & ~isempty(dof.(dir1)) )
            u.(vel1)(dof.(dir1)) = dirichlet.(ddir1);
        end
        if ( isfield(dof,dir2) & ~isempty(dof.(dir2)) )
            u.(vel2)(dof.(dir2)) = dirichlet.(ddir2);
        end
    end
elseif settings.darcymixed>3
    % Essential boundary conditions are on the pressure
    for i=1:length(domain)
        dname         = strcat('d',domain(i));
        internal = strcat(dname,'_internal');
        dir1    = strcat(dname,'_dirichlet_3');
        ddir1   = strcat(dname,'_3');
        nunodes = strcat(dname,'_u');
        npnodes = strcat(dname,'_p');
        vel1    = strcat(dname,'_1');
        vel2    = strcat(dname,'_2');
        press   = strcat(dname,'_3');
        % Allocate memory space
        u.(vel1)  = zeros(2*dof.(nunodes),1);
        u.(vel2)  = zeros(dof.(nunodes),1);
        u.(press) = zeros(dof.(npnodes),1);
        %
        % Velocity
        u.(vel1) = un(1:dof.(nunodes));
        u.(vel2) = un(dof.(nunodes)+1:2*dof.(nunodes));
        % Pressure
        u.(press)(dof.(internal)(2:end),1) = un(2*dof.(nunodes)+1:end);
        clear un;
        % Set Dirichlet values
        if ( isfield(dof,dir1) & ~isempty(dof.(dir1)) )
            u.(press)(dof.(dir1)) = dirichlet.(ddir1);
        end
    end
end

return