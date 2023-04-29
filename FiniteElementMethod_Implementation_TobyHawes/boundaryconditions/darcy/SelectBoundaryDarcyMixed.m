function [dof,matdata] = SelectBoundaryDarcyMixed (mesh,dof,bd_data,...
    domain,matdata,settings)

% Select type of boundary for the Darcy problem in mixed form.

bd_flag = ['1' '2']; % first or second component of the velocity
                     % The second component is selected only for mixed bc
                     % on the same edge

dname          = strcat('d',domain);
edges_name     = strcat(dname,'_e');
dirichlet_name = strcat(dname,'_DirichletEdges');
neumann_name   = strcat(dname,'_NeumannEdges');
nodes_name_v   = strcat(dname,'_all_u');
nodes_name_p   = strcat(dname,'_all_p');
unodes         = strcat(dname,'_u');
fem            = strcat(dname,'_fem');
bctype         = strcat(dname,'_bctype');

switch mesh.(fem)
    case 'P11'
        vindk = 2;
        pindk = 2;
        inde  = 5;
        nx    = 3;
    case 'P12'
        vindk = 2;
        pindk = 3;
        inde  = 6;
        nx    = 4;
    case 'P21'
        vindk = 3;
        pindk = 2;
        inde  = 6;
        nx    = 4;
    case 'P22'
        vindk = 3;
        pindk = 3;
        inde  = 6;
        nx    = 4;
end

% This flag takes into account whether we should fix the value of the
% pressure in one point
flag = 0;
dof.(bctype) = [0 0 0];

if settings.darcymixed<4
    % Normal velocity is the essential variable
    % Pressure is the natural variable
    % Search Dirichlet (essential) nodes
    if ~isempty(bd_data.(dirichlet_name))
        dof.(bctype)(1) = 1;
        edges = zeros(1,size(bd_data.(dirichlet_name),2));
        m = 0;
        for k=1:size(bd_data.(dirichlet_name),2)
            if bd_data.(dirichlet_name)(1,k)~=0
                m = m+1;
                edges(m) = bd_data.(dirichlet_name)(1,k);
            end
        end
        edges = edges(1,1:m);
        indices1 = zeros(1,size(mesh.(edges_name),2));
        indices2 = indices1;
        m1 = 0;
        m2 = 0;
        for j=1:length(edges)
            for k=1:size(mesh.(edges_name),2)
                if mesh.(edges_name)(inde,k)==edges(j)
                    if mesh.(edges_name)(nx,k)~=0
                        m1 = m1 + 1;
                        indices1(1,m1) = k;
                    else
                        m2 = m2 + 1;
                        indices2(1,m2) = k;
                    end
                end
            end
        end
        indices1 = indices1(1,1:m1);
        indices2 = indices2(1,1:m2);
        % Dirichlet nodes (1st velocity component)
        diri1 = strcat(dname,'_dirichlet_1');
        inti1 = strcat(dname,'_internal_1');
        dof.(diri1) = mesh.(edges_name)(1:vindk,indices1);
        for j=1:size(dof.(diri1),2)
            dof.(diri1)(:,j) = (mesh.(edges_name)(nx,indices1(1,j)))*dof.(diri1)(:,j);
        end
        dof.(diri1) = reshape(dof.(diri1),1,vindk*m1);
        dof.(diri1) = unique(dof.(diri1));
        dof.(inti1) = setdiff(dof.(nodes_name_v),abs(dof.(diri1)));
        % Dirichlet nodes (2nd velocity component)
        diri2 = strcat(dname,'_dirichlet_2');
        inti2 = strcat(dname,'_internal_2');
        dof.(diri2) = mesh.(edges_name)(1:vindk,indices2);
        for j=1:size(dof.(diri2),2)
            dof.(diri2)(:,j) = (mesh.(edges_name)(nx+1,indices2(1,j)))*dof.(diri2)(:,j);
        end
        dof.(diri2) = reshape(dof.(diri2),1,vindk*m2);
        dof.(diri2) = unique(dof.(diri2));
        dof.(inti2) = setdiff(dof.(nodes_name_v),abs(dof.(diri2)));
        %
    else
        for i=1:length(bd_flag)
            diri = strcat(dname,'_dirichlet_',bd_flag(i));
            inti = strcat(dname,'_internal_',bd_flag(i));
            dof.(diri) = [];
            dof.(inti) = dof.(nodes_name_v);
        end
    end
    % Compact internal nodes
    inti1 = strcat(dname,'_internal_1');
    inti2 = strcat(dname,'_internal_2');
    intflag = 0;
    if ~isfield(dof,inti2)
        inti2 = inti1;
        intflag = 1;
    end
    inti  = strcat(dname,'_internal');
    nint1 = length(dof.(inti1));
    nint2 = length(dof.(inti2));
    dof.(inti) = [nint1+2, nint2+nint1+2, dof.(inti1), dof.(inti2)+dof.(unodes)];
    dof = rmfield(dof,inti1);
    if ~intflag
        dof = rmfield(dof,inti2);
    end
    % Search Neumann (natural) boundaries
    if ~isempty(bd_data.(neumann_name))
        dof.(bctype)(2) = 1;
        neui = strcat(dname,'_neumann_1');
        edges = zeros(1,size(bd_data.(neumann_name),2));
        m = 0;
        for k=1:size(bd_data.(neumann_name),2)
            if bd_data.(neumann_name)(1,k)~=0
                m = m+1;
                edges(m) = bd_data.(neumann_name)(1,k);
            end
        end
        edges = edges(1,1:m);
        indices = zeros(1,size(mesh.(edges_name),2));
        m = 0;
        for j=1:length(edges)
            for k=1:size(mesh.(edges_name),2)
                if mesh.(edges_name)(inde,k)==edges(j)
                    m = m+1;
                    indices(m) = k;
                end
            end
        end
        indices = indices(1,1:m);
        dof.(neui) = unique(indices);
        %
        flag = 1;
    else
        neui = strcat(dname,'_neumann_1');
        dof.(neui) = [];
    end
    %
elseif settings.darcymixed>3
    % Normal velocity is the natural variable
    % Pressure is the essential variable
    % Search Dirichlet (essential) nodes [pressure]
    if ~isempty(bd_data.(dirichlet_name))
        dof.(bctype)(1) = 1;
        edges = zeros(1,size(bd_data.(dirichlet_name),2));
        m = 0;
        for k=1:size(bd_data.(dirichlet_name),2)
            if bd_data.(dirichlet_name)(1,k)~=0
                m = m+1;
                edges(m) = bd_data.(dirichlet_name)(1,k);
            end
        end
        edges = edges(1,1:m);
        indices1 = zeros(1,size(mesh.(edges_name),2));
        m1 = 0;
        for j=1:length(edges)
            for k=1:size(mesh.(edges_name),2)
                if mesh.(edges_name)(inde,k)==edges(j)
                    m1 = m1 + 1;
                    indices1(1,m1) = k;
                end
            end
        end
        indices1 = indices1(1,1:m1);
        % Dirichlet nodes (pressure)
        diri1 = strcat(dname,'_dirichlet_3');
        inti1 = strcat(dname,'_internal_3');
        dof.(diri1) = mesh.(edges_name)(1:pindk,indices1);
        dof.(diri1) = reshape(dof.(diri1),1,pindk*m1);
        dof.(diri1) = unique(dof.(diri1));
        dof.(inti1) = setdiff(dof.(nodes_name_p),dof.(diri1));
        %
    else
        diri = strcat(dname,'_dirichlet_3');
        inti = strcat(dname,'_internal_3');
        dof.(diri) = [];
        dof.(inti) = dof.(nodes_name_p);
    end
    % Compact internal nodes
    inti1 = strcat(dname,'_internal_3');
    inti  = strcat(dname,'_internal');
    nint1 = length(dof.(inti1));
    dof.(inti) = [nint1+1, dof.(inti1)];
    dof = rmfield(dof,inti1);
    % Search Neumann (natural) boundaries [normal velocity]
    if ~isempty(bd_data.(neumann_name))
        dof.(bctype)(2) = 1;
        neui = strcat(dname,'_neumann_1');
        edges = zeros(1,size(bd_data.(neumann_name),2));
        m = 0;
        for k=1:size(bd_data.(neumann_name),2)
            if bd_data.(neumann_name)(1,k)~=0
                m = m+1;
                edges(m) = bd_data.(neumann_name)(1,k);
            end
        end
        edges = edges(1,1:m);
        indices = zeros(1,size(mesh.(edges_name),2));
        m = 0;
        for j=1:length(edges)
            for k=1:size(mesh.(edges_name),2)
                if mesh.(edges_name)(inde,k)==edges(j)
                    m = m+1;
                    indices(m) = k;
                end
            end
        end
        indices = indices(1,1:m);
        dof.(neui) = unique(indices);
        %
        flag = 1;
    else
        neui = strcat(dname,'_neumann_1');
        dof.(neui) = [];
    end
    %
end

if ~flag
    %
    fprintf('   [SelectDarcyMixedStokes] Pressure defined up to an additive constant\n');
    fprintf('                            (null average is automatically imposed)\n');
    %
end

return