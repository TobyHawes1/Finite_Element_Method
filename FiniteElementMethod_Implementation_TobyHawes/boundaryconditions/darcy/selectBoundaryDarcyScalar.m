function [dof] = selectBoundaryDarcyScalar (dof,mesh,boundaryData,domain)

bdFlag = {'3'}; % The unknown is scalar (pressure)

dName         = strcat('d',domain);
edgesName     = strcat(dName,'_e');
dirichletName = strcat(dName,'_DirichletEdges');
neumannName   = strcat(dName,'_NeumannEdges');
robinName     = strcat(dName,'_RobinEdges');
nodesName     = strcat(dName,'_all_v3');
fem           = strcat(dName,'_fem');
bcType        = strcat(dName,'_bctype');

switch mesh.(fem)
    case 'P1P1'
        indk = 2;
        inde = 5;
    case 'P1P2'
        indk = 3;
        inde = 6;
    case 'P2P1'
        indk = 2;
        inde = 6;
    case 'P2P2'
        indk = 3;
        inde = 6;
    case 'Q1Q1'
        indk = 2;
        inde = 5;
    case 'Q1Q2'
        indk = 3;
        inde = 6;
    case 'Q2Q1'
        indk = 2;
        inde = 6;
    case 'Q2Q2'
        indk = 3;
        inde = 6;
end

dof.(bcType) = [0 0 0];

% Searching Dirichlet nodes
diri = strcat(dName,'_dirichlet_v',bdFlag{1});
inti = strcat(dName,'_internal');
nint = strcat(dName,'_ninternal');
dof.(diri) = [];
if ~isempty(boundaryData.(dirichletName))
    dof.(bcType)(1) = 1;
    edges = zeros(1,size(boundaryData.(dirichletName),2));
    m = 0;
    for k=1:size(boundaryData.(dirichletName),2)
        if boundaryData.(dirichletName)(1,k)~=0
            m = m+1;
            edges(m) = boundaryData.(dirichletName)(1,k);
        end
    end
    edges = edges(1,1:m);
    indices = zeros(1,size(mesh.(edgesName),2));
    m = 0;
    for j=1:length(edges)
        for k=1:size(mesh.(edgesName),2)
            if mesh.(edgesName)(inde,k)==edges(j)
                m = m+1;
                indices(m) = k;
            end
        end
    end
    indices = indices(1,1:m);
    dof.(diri) = reshape(mesh.(edgesName)(1:indk,indices),1,indk*m);
    dof.(diri) = unique(dof.(diri));
    dof.(inti) = setdiff(dof.(nodesName),dof.(diri));
    dof.(nint) = length(dof.(inti));
else
    dof.(inti) = dof.(nodesName);
    dof.(nint) = length(dof.(inti));
end
dof.(inti) = [dof.(nint)+1, dof.(inti)];
dof = rmfield(dof,nint);


% Searching Neumann boundaries
neui = strcat(dName,'_neumann_v',bdFlag{1});
if ~isempty(boundaryData.(neumannName))
    dof.(bcType)(2) = 1;
    edges = zeros(1,size(boundaryData.(neumannName),2));
    m = 0;
    for k=1:size(boundaryData.(neumannName),2)
        if boundaryData.(neumannName)(1,k)~=0
            m = m+1;
            edges(m) = boundaryData.(neumannName)(1,k);
        end
    end
    edges = edges(1,1:m);
    indices = zeros(1,size(mesh.(edgesName),2));
    m = 0;
    for j=1:length(edges)
        for k=1:size(mesh.(edgesName),2)
            if mesh.(edgesName)(inde,k)==edges(j)
                m = m+1;
                indices(m) = k;
            end
        end
    end
    indices = indices(1,1:m);
    dof.(neui) = unique(indices);
else
    dof.(neui) = [];
end

% Searching Robin boundaries
robi = strcat(dName,'_robin_v',bdFlag{1});
if ~isempty(boundaryData.(robinName))
    dof.(bcType)(3) = 1;
    edges = zeros(1,size(boundaryData.(robinName),2));
    m = 0;
    for k=1:size(boundaryData.(robinName),2)
        if boundaryData.(robinName)(1,k)~=0
            m = m+1;
            edges(m) = boundaryData.(robinName)(1,k);
        end
    end
    edges = edges(1,1:m);
    indices = zeros(1,size(mesh.(edgesName),2));
    m = 0;
    for j=1:length(edges)
        for k=1:size(mesh.(edgesName),2)
            if mesh.(edgesName)(inde,k)==edges(j)
                m = m+1;
                indices(m) = k;
            end
        end
    end
    indices = indices(1,1:m);
    dof.(robi) = unique(indices);
else
    dof.(robi) = [];
end

end