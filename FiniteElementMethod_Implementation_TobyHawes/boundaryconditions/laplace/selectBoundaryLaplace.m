function [dof] = selectBoundaryLaplace (dof,mesh,boundaryData,domain,variable)

% bdFlag = {'v1'}; % the unknown is scalar (one component)

% If the domain is the coarse one, then only Dirichlet nodes are selected
% (homogeneous Neumann conditions are imposed elsewhere, since they will be
% used in the framework of preconditioning)

domainVariable = strcat(domain,'_',variable);
edgesName      = strcat(domainVariable,'_e');
dirichletName  = strcat(domain,'_DirichletEdges');
neumannName    = strcat(domain,'_NeumannEdges');
robinName      = strcat(domain,'_RobinEdges');
fem            = strcat(domainVariable,'_fem');
bcType         = strcat(domainVariable,'_bctype');
fem = mesh.(fem);
%fem = fem(end-1:end);

% Set type of elements
switch fem
    case {'Q1','P1'}
        nodesOnEdge = 2;
    case {'Q2','P2'}
        nodesOnEdge = 3;
    case {'Q3','Q3gl','P3','P3gl'}
        nodesOnEdge = 4;
end

dof.(bcType) = [0 0 0];

% Searching Dirichlet nodes
diri = strcat(domain,'_',variable,'_1_dirichlet');
inti = strcat(domain,'_',variable,'_1_internal');
nint = strcat(domain,'_',variable,'_1_ninternal');
dof.(diri) = [];
if ~isempty(boundaryData.(dirichletName))
    dof.(bcType)(1) = 1;
    %    for i = 1:size(bdData.(dirichletName),1)
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
            if mesh.(edgesName)(1,k)==edges(j)
                m = m+1;
                indices(m) = k;
            end
        end
    end
    indices = indices(1,1:m);
    dof.(diri) = reshape(mesh.(edgesName)(4:3+nodesOnEdge,indices),1,nodesOnEdge*m);
    dof.(diri) = unique(dof.(diri));
    dof.(inti) = setdiff(dof.(domainVariable),dof.(diri));
    dof.(nint) = length(dof.(inti));
    %    end
else
    dof.(inti) = dof.(domainVariable);
    dof.(nint) = length(dof.(inti));
end

% Searching Neumann boundaries
neui = strcat(domain,'_',variable,'_1_neumann');
if (~isempty(boundaryData.(neumannName)) && ~strcmp(domain,'C'))
    dof.(bcType)(2) = 1;
    %     for i=1:size(bdData.(neumannName),1)
    %         neui = strcat(dName,'_neumann_',bdFlag{i});
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
            if mesh.(edgesName)(1,k)==edges(j)
                m = m+1;
                indices(m) = k;
            end
        end
    end
    indices = indices(1,1:m);
    dof.(neui) = unique(indices);
    %     end
else
    %     for i=1:length(bdFlag)
    %         neui = strcat(dName,'_neumann_',bdFlag{i});
    dof.(neui) = [];
    %     end
end

% Searching Robin boundaries
robi = strcat(domain,'_',variable,'_1_robin');
if (~isempty(boundaryData.(robinName)) && ~strcmp(domain,'C'))
    dof.(bcType)(3) = 1;
    %     for i=1:size(bdData.(robinName),1)
    %         robi = strcat(dName,'_robin_',bdFlag{i});
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
            if mesh.(edgesName)(1,k)==edges(j)
                m = m+1;
                indices(m) = k;
            end
        end
    end
    indices = indices(1,1:m);
    dof.(robi) = unique(indices);
    %     end
else
    %     for i=1:length(bdFlag)
    %         robi = strcat(dName,'_robin_',bdFlag{i});
    dof.(robi) = [];
    %     end
end

end