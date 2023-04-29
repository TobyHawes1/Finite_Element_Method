function [dof] = selectParametricBoundaryStokes (mesh,dof,bdData,domain)

% Select type of boundary for Stokes problem.

velocityComponents = 2;

domainVelocity = strcat(domain,'_v1');
edgesVelocity  = strcat(domainVelocity,'_e');
dirichletName  = strcat(domain,'_ParametricDirichletEdges');
neumannName    = strcat(domain,'_ParametricNeumannEdges');
robinName      = strcat(domain,'_ParametricRobinEdges');
femVelocity    = strcat(domainVelocity,'_fem');
femVelocity    = mesh.(femVelocity);
bcType         = strcat(domain,'_parametricbctype');

switch femVelocity
    case {'Q1'}%'P1'
        nodesOnEdge = 2;
    case {'Q2'}%'P2'
        nodesOnEdge = 3;
    case {'Q3','Q3gl'}
        nodesOnEdge = 4;
end

dof.(bcType) = [0 0 0];

% Searching Dirichlet nodes
if ~isempty(bdData.(dirichletName))
    dof.(bcType)(1) = 1;
    for i = 1:size(bdData.(dirichletName),1)
        dirichlet = strcat(domainVelocity,'_',int2str(i),'_parametricDirichlet');
        dof.(dirichlet) = [];
        edges = zeros(1,size(bdData.(dirichletName),2));
        m = 0;
        for k = 1:size(bdData.(dirichletName),2)
            if bdData.(dirichletName)(i,k)~=0
                m = m+1;
                edges(m) = bdData.(dirichletName)(i,k);
            end
        end
        edges = edges(1,1:m);
        indices = zeros(1,size(mesh.(edgesVelocity),2));
        m = 0;
        for j = 1:length(edges)
            for k = 1:size(mesh.(edgesVelocity),2)
                if mesh.(edgesVelocity)(1,k) == edges(j)
                    m = m+1;
                    indices(m) = k;
                end
            end
        end
        indices = indices(1,1:m);
        dof.(dirichlet) = reshape(mesh.(edgesVelocity)(4:3+nodesOnEdge,indices),1,nodesOnEdge*m);
        dof.(dirichlet) = unique(dof.(dirichlet));
    end
    %
    if size(bdData.(dirichletName),1) == 1
        dirichlet1 = strcat(domainVelocity,'_1_parametricDirichlet');
        dirichlet2 = strcat(domainVelocity,'_2_parametricDirichlet');
        dof.(dirichlet2) = dof.(dirichlet1);
    end
    %
else
    for i = 1:velocityComponents
        dirichlet = strcat(domainVelocity,'_',int2str(i),'parametricDirichlet');
        dof.(dirichlet) = [];
    end
end

% Searching Neumann boundaries
if ~isempty(bdData.(neumannName))
    neumann = strcat(domainVelocity,'_parametricNeumann');
    edges = zeros(1,size(bdData.(neumannName),2));
    m = 0;
    for k=1:size(bdData.(neumannName),2)
        if bdData.(neumannName)(1,k)~=0
            m = m+1;
            edges(m) = bdData.(neumannName)(1,k);
        end
    end
    dof.(bcType)(2) = m;
    edges = edges(1,1:m);
    for j=1:length(edges)
        indices = zeros(1,size(mesh.(edgesVelocity),2));
        m = 0;
        for k=1:size(mesh.(edgesVelocity),2)
            if mesh.(edgesVelocity)(1,k)==edges(j)
                m = m+1;
                indices(m) = k;
            end
        end
        indices = indices(1,1:m);
        neumannEdge = strcat(neumann,'_edge_',int2str(j)); %(edges(1,j)));
        dof.(neumannEdge) = unique(indices);
    end
else
    neumann = strcat(domainVelocity,'_parametricNeumann');
    dof.(neumann) = [];
end

% Searching Robin boundaries
%if (~isempty(bdData.(robinName)) && ~strcmp(domain,'C'))
if ~isempty(bdData.(robinName))
    dof.(bcType)(3) = 1;
    for i=1:size(bdData.(robinName),1)
        robin = strcat(domainVelocity,'_',int2str(i),'_parametricRobin');
        edges = zeros(1,size(bdData.(robinName),2));
        m = 0;
        for k=1:size(bdData.(robinName),2)
            if bdData.(robinName)(i,k)~=0
                m = m+1;
                edges(m) = bdData.(robinName)(i,k);
            end
        end
        edges = edges(1,1:m);
        indices = zeros(1,size(mesh.(edgesVelocity),2));
        m = 0;
        for j=1:length(edges)
            for k=1:size(mesh.(edgesVelocity),2)
                if mesh.(edgesVelocity)(1,k)==edges(j)
                    m = m+1;
                    indices(m) = k;
                end
            end
        end
        indices = indices(1,1:m);
        dof.(robin) = unique(indices);
    end
    flag = 1;
    %
    if size(bdData.(robinName),1) == 1
        robin1 = strcat(domainVelocity,'_1_parametricRobin');
        robin2 = strcat(domainVelocity,'_2_parametricRobin');
        dof.(robin2) = dof.(robin1);
    end
    %
else
    for i = 1:velocityComponents
        robin = strcat(domainVelocity,'_',int2str(i),'_parametricRobin');
        dof.(robin) = [];
    end
end

end