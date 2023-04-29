function [dof] = selectBoundaryStokes (mesh,dof,bdData,domain)

% Select type of boundary for Stokes problem.

velocityComponents = 2;
% If the domain is the coarse one, then only Dirichlet nodes are selected
% (homogeneous Neumann conditions are imposed elsewhere, since they will be
% used in the framework of preconditioning)

domainVelocity = strcat(domain,'_v1');
edgesVelocity  = strcat(domainVelocity,'_e');
dirichletName  = strcat(domain,'_DirichletEdges');
neumannName    = strcat(domain,'_NeumannEdges');
robinName      = strcat(domain,'_RobinEdges');
numberOfNodesVelocity = strcat(domain,'_nv1');
femVelocity    = strcat(domainVelocity,'_fem');
femVelocity    = mesh.(femVelocity);
bcType         = strcat(domain,'_bctype');

switch femVelocity
    case {'Q1'}%'P1'
        nodesOnEdge = 2;
    case {'Q2'}%'P2'
        nodesOnEdge = 3;
    case {'Q3','Q3gl'}
        nodesOnEdge = 4;
end

% This flag takes into account whether we should fix the value of the
% pressure in one point
flag = 0;
dof.(bcType) = [0 0 0];

% Searching Dirichlet nodes
if ~isempty(bdData.(dirichletName))
    dof.(bcType)(1) = 1;
    for i = 1:size(bdData.(dirichletName),1)
        dirichlet = strcat(domainVelocity,'_',int2str(i),'_dirichlet');
        internal = strcat(domainVelocity,'_',int2str(i),'_internal');
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
        dof.(internal) = setdiff(dof.(domainVelocity),dof.(dirichlet));
    end
    %
    if size(bdData.(dirichletName),1) == 1
        dirichlet1 = strcat(domainVelocity,'_1_dirichlet');
        internal1 = strcat(domainVelocity,'_1_internal');
        dirichlet2 = strcat(domainVelocity,'_2_dirichlet');
        internal2 = strcat(domainVelocity,'_2_internal');
        dof.(dirichlet2) = dof.(dirichlet1);
        dof.(internal2) = dof.(internal1);
    end
    %
else
    for i = 1:velocityComponents
        dirichlet = strcat(domainVelocity,'_',int2str(i),'_dirichlet');
        internal = strcat(domainVelocity,'_',int2str(i),'_internal');
        dof.(dirichlet) = [];
        dof.(internal) = dof.(domainVelocity);
    end
end
% % Gather information on internal nodes
internal1  = strcat(domainVelocity,'_1_internal');
internal2  = strcat(domainVelocity,'_2_internal');
internal   = strcat(domainVelocity,'_internal');
dof.(internal) = [dof.(internal1), dof.(internal2)+dof.(numberOfNodesVelocity)];

numberOfInternal1 = strcat(domainVelocity,'_1_ninternal');
numberOfInternal2 = strcat(domainVelocity,'_2_ninternal');
numberOfInternal  = strcat(domainVelocity,'_ninternal');

dof.(numberOfInternal1) = size(dof.(internal1),2);
dof.(numberOfInternal2) = size(dof.(internal2),2);
dof.(numberOfInternal) = dof.(numberOfInternal1) + dof.(numberOfInternal2);

% Searching Neumann boundaries
%if (~isempty(bdData.(neumannName)) && ~strcmp(domain,'C'))
if ~isempty(bdData.(neumannName))
    dof.(bcType)(2) = 1;
    for i = 1:size(bdData.(neumannName),1)
        neumann = strcat(domainVelocity,'_',int2str(i),'_neumann');
        edges = zeros(1,size(bdData.(neumannName),2));
        m = 0;
        for k=1:size(bdData.(neumannName),2)
            if bdData.(neumannName)(i,k)~=0
                m = m+1;
                edges(m) = bdData.(neumannName)(i,k);
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
        dof.(neumann) = unique(indices);
    end
    flag = 1;
    %
    if size(bdData.(neumannName),1) == 1
        neumann1 = strcat(domainVelocity,'_1_neumann');
        neumann2 = strcat(domainVelocity,'_2_neumann');
        dof.(neumann2) = dof.(neumann1);
    end
    %
else
    for i = 1:velocityComponents
        neumann = strcat(domainVelocity,'_',int2str(i),'_neumann');
        dof.(neumann) = [];
    end
end

% Searching Robin boundaries
%if (~isempty(bdData.(robinName)) && ~strcmp(domain,'C'))
if ~isempty(bdData.(robinName))
    dof.(bcType)(3) = 1;
    for i=1:size(bdData.(robinName),1)
        robin = strcat(domainVelocity,'_',int2str(i),'_robin');
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
        robin1 = strcat(domainVelocity,'_1_robin');
        robin2 = strcat(domainVelocity,'_2_robin');
        dof.(robin2) = dof.(robin1);
    end
    %
else
    for i = 1:velocityComponents
        robin = strcat(domainVelocity,'_',int2str(i),'_robin');
        dof.(robin) = [];
    end
end

if ~flag
    fprintf('   [SelectBoundaryStokes] Pressure defined up to an additive constant\n');
    fprintf('                          (null average is automatically imposed)\n');
end

end