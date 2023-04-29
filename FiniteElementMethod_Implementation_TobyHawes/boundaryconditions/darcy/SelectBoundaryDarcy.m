function [dof] = SelectBoundaryDarcy (dof,mesh,bd_data,domain)

bd_flag = '3'; % the unknown is scalar (pressure)

dname          = strcat('d',domain);
edges_name     = strcat(dname,'_e');
dirichlet_name = strcat(dname,'_Dboundary');
neumann_name   = strcat(dname,'_Nboundary');
robin_name     = strcat(dname,'_Rboundary');
nodes_name     = strcat(dname,'_all_v3');
fem            = strcat(dname,'_fem');
bctype         = strcat(dname,'_bctype');

% Set type of elements
if ( strcmp(mesh.(fem),'P22') | strcmp(mesh.(fem),'P12') )
    indk = 3;
    inde = 6;
elseif strcmp(mesh.(fem),'P21')
    indk = 2;
    inde = 6;
elseif strcmp(mesh.(fem),'P11')
    indk = 2;
    inde = 5;
end

dof.(bctype) = [0 0 0];

% Searching Dirichlet nodes
if ~isempty(bd_data.(dirichlet_name))
    dof.(bctype)(1) = 1;
    for i=1:size(bd_data.(dirichlet_name),1)
        diri = strcat(dname,'_dirichlet_',bd_flag(i));
        inti = strcat(dname,'_internal_',bd_flag(i));
        nint = strcat(dname,'_ninternal_',bd_flag(i));
        dof.(diri) = [];
        edges = zeros(1,size(bd_data.(dirichlet_name),2));
        m = 0;
        for k=1:size(bd_data.(dirichlet_name),2)
            if bd_data.(dirichlet_name)(i,k)~=0
                m = m+1;
                edges(m) = bd_data.(dirichlet_name)(i,k);
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
        dof.(diri) = reshape(mesh.(edges_name)(1:indk,indices),1,indk*m);
        dof.(diri) = unique(dof.(diri));
        dof.(inti) = setdiff(dof.(nodes_name),dof.(diri));
        dof.(nint) = length(dof.(inti));
    end
else
    diri = strcat(dname,'_dirichlet_',bd_flag(1));
    inti = strcat(dname,'_internal_',bd_flag(1));
    nint = strcat(dname,'_ninternal_',bd_flag(1));
    dof.(diri) = [];
    dof.(inti) = dof.(nodes_name);
    dof.(nint) = length(dof.(inti));
end

% Searching Neumann boundaries
if (~isempty(bd_data.(neumann_name)) & ~strcmp(domain,'C'))
    dof.(bctype)(2) = 1;
    for i=1:size(bd_data.(neumann_name),1)
        neui = strcat(dname,'_neumann_',bd_flag(i));
        edges = zeros(1,size(bd_data.(neumann_name),2));
        m = 0;
        for k=1:size(bd_data.(neumann_name),2)
            if bd_data.(neumann_name)(i,k)~=0
                m = m+1;
                edges(m) = bd_data.(neumann_name)(i,k);
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
    end
else
    for i=1:length(bd_flag)
        neui = strcat(dname,'_neumann_',bd_flag(i));
        dof.(neui) = [];
    end
end

% Searching Robin boundaries
if (~isempty(bd_data.(robin_name)) & ~strcmp(domain,'C'))
    dof.(bctype)(3) = 1;
    for i=1:size(bd_data.(robin_name),1)
        robi = strcat(dname,'_robin_',bd_flag(i));
        edges = zeros(1,size(bd_data.(robin_name),2));
        m = 0;
        for k=1:size(bd_data.(robin_name),2)
            if bd_data.(robin_name)(i,k)~=0
                m = m+1;
                edges(m) = bd_data.(robin_name)(i,k);
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
        dof.(robi) = unique(indices);
    end
else
    for i=1:length(bd_flag)
        robi = strcat(dname,'_robin_',bd_flag(i));
        dof.(robi) = [];
    end
end

return