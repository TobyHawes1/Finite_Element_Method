function [rhs,varargout] = stokesNaturalBC (rhs,mesh,dof,domain,data,boundaryType,time,varargin)

% StokesNaturalBC imposes natural boundary conditions to the Stokes
%    equations. The boundary conditions can be of Neumann or Robin type.
%
% 'type' is a string which can take the values 'neumann' or 'robin'.
%
% In case of NEUMANN conditions, we need to define the two components of
% (\nu\nabla\mathbf{u}\cdot\mathbf{n}-p\mathbf{n})
%
% In case of ROBIN conditions, we need to define the two components of
% (\nu\nabla\mathbf{u}\cdot\mathbf{n}-p\mathbf{n}) + a * \mathbf{u}

component = zeros(1,2); m = 0;
domainVelocity = strcat(domain,'_v1');
nodesVelocity  = strcat(domainVelocity,'_p');
edgesVelocity  = strcat(domainVelocity,'_e');
numberOfNodesVelocity = strcat(domain,'_nv1');
femVelocity    = strcat(domainVelocity,'_fem');
femVelocity    = mesh.(femVelocity);

if strcmp(boundaryType,'neumann')
    dataFlag = 22; % Neumann
elseif strcmp(boundaryType,'robin')
    dataFlag = 23; % Robin
end

if ~isempty(varargin)
    arguments = varargin{1};
else
    arguments = [];
end

% Load quadrature nodes and weights
[ti,wi,phi] = basisOnQuad1D(femVelocity,10);
switch femVelocity
    case 'Q1'
        nodesOnEdge = 2;
    case 'Q2'
        nodesOnEdge = 3;
    case {'Q3','Q3gl'}
        nodesOnEdge = 4;
end

% First component of the velocity
naturalName = strcat(domainVelocity,'_1_',boundaryType);
for i = 1:size(dof.(naturalName),2)
    x = mesh.(nodesVelocity)(1,mesh.(edgesVelocity)(4:5,dof.(naturalName)(i)));
    y = mesh.(nodesVelocity)(2,mesh.(edgesVelocity)(4:5,dof.(naturalName)(i)));
    edgeLength = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    normalVector = mesh.(edgesVelocity)([2;3],dof.(naturalName)(i));
    % Transformation to the reference interval [0,1]
    %x = ((1-ti)*x(1)+(ti)*x(2));
    %y = ((1-ti)*y(1)+(ti)*y(2));
    %fnxy = feval(data,dataFlag,x,y,time,1,normalVector,arguments);
    %ind = mesh.(edgesVelocity)(4:3+nodesOnEdge,dof.(naturalName)(i));
    %rhs.(domainVelocity)(ind) = rhs.(domainVelocity)(ind) + edgeLength*((wi.*fnxy)*phi')';
    % Transformation to the reference interval [-1,1]
    x = (x(2)-x(1))*ti/2 + (x(2)+x(1))/2;
    y = (y(2)-y(1))*ti/2 + (y(2)+y(1))/2;
    fnxy = feval(data,dataFlag,x,y,time,1,normalVector,arguments);
    ind = mesh.(edgesVelocity)(4:3+nodesOnEdge,dof.(naturalName)(i));
    rhs.(domainVelocity)(ind) = rhs.(domainVelocity)(ind) + ...
        (edgeLength/2)*((wi.*fnxy)*phi')';
    if ( strcmp(boundaryType,'robin') && i==1 )
        m = m+1;
        component(m) = 1;
    end
end

% Second component of the velocity
naturalName = strcat(domainVelocity,'_2_',boundaryType);
% if ~isfield(dof,naturalName)
%     naturalName = strcat(dName,'_',boundaryType,'_v1');
% end
for i = 1:size(dof.(naturalName),2)
    x = mesh.(nodesVelocity)(1,mesh.(edgesVelocity)(4:5,dof.(naturalName)(i)));
    y = mesh.(nodesVelocity)(2,mesh.(edgesVelocity)(4:5,dof.(naturalName)(i)));
    edgeLength = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    normalVector = mesh.(edgesVelocity)([2;3],dof.(naturalName)(i));
    % Transformation to the reference interval [0,1]
    %x = ((1-ti)*x(1)+(ti)*x(2));
    %y = ((1-ti)*y(1)+(ti)*y(2));
    % Transformation to the reference interval [-1,1]
    x = (x(2)-x(1))*ti/2 + (x(2)+x(1))/2;
    y = (y(2)-y(1))*ti/2 + (y(2)+y(1))/2;
    fnxy = feval(data,dataFlag,x,y,time,2,normalVector,arguments);
    ind = mesh.(edgesVelocity)(4:3+nodesOnEdge,dof.(naturalName)(i));
    rhs.(domainVelocity)(ind+dof.(numberOfNodesVelocity)) = ...
        rhs.(domainVelocity)(ind+dof.(numberOfNodesVelocity)) + ...
        (edgeLength/2)*((wi.*fnxy)*phi')';
    if ( strcmp(boundaryType,'robin') && i==1 )
        m = m+1;
        component(m) = 2;
    end
end
if m>0
    varargout(1) = {component(1:m)};
else
    varargout(1) = {[]};
end
    
end