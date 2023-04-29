function [rhs] = laplaceNaturalBC (rhs,mesh,dof,domain,variable,data,type,time,varargin)

% LaplaceNaturalBC imposes natural boundary conditions to the Poisson
%    equations. The boundary conditions can be of Neumann or Robin type.
%
% 'type' is a string which can take the values 'neumann' or 'robin'.
%
% In case of NEUMANN conditions, we need to define the conormal derivative
% \nu \nabla u \cdot \mathbf{n}
%
% In case of ROBIN conditions, we need to define the function
% \nu \nabla u \cdot \mathbf{n} + a * u

domainVariable = strcat(domain,'_',variable);
edges    = strcat(domainVariable,'_e');
nodes    = strcat(domainVariable,'_p');
fem      = strcat(domainVariable,'_fem');
fem = mesh.(fem);
%fem = fem(end-1:end);
naturalName = strcat(domainVariable,'_1_',type);

% Load quadrature nodes and weights
[ti,wi,phi] = basisOnQuad1D(fem,10);
switch fem
    case {'P1','Q1'}
        nodesOnEdge = 2;
    case {'P2','Q2'}
        nodesOnEdge = 3;
    case {'P3','P3gl','Q3','Q3gl'}
        nodesOnEdge = 4;
end

if strcmp(type,'neumann')
    dataFlagNatural = 22; % Neumann
elseif strcmp(type,'robin')
    dataFlagNatural = 23; % Robin
end

for i=1:length(dof.(naturalName))
    x = mesh.(nodes)(1,mesh.(edges)(4:5,dof.(naturalName)(i)));
    y = mesh.(nodes)(2,mesh.(edges)(4:5,dof.(naturalName)(i)));
    edgeLength = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    normalVector = mesh.(edges)([2;3],dof.(naturalName)(i));
    % Transformation to the reference interval [0,1]
    % x = ((1-ti)*x(1)+(ti)*x(2));
    % y = ((1-ti)*y(1)+(ti)*y(2));
    % fnxy = feval(data,dataFlag,x,y,normalVector,varargin{:});
    % ind = mesh.(edges)(4:3+nodesOnEdge,dof.(naturalName)(i));
    % rhs.(domainVariable)(ind) = rhs.(domainVariable)(ind) + ...
    %     edgeLength*((wi.*fnxy)*phi')';
    % Transformation to the reference interval [-1,1]
    x = (x(2)-x(1))*ti/2 + (x(2)+x(1))/2;
    y = (y(2)-y(1))*ti/2 + (y(2)+y(1))/2;
    fnxy = feval(data,dataFlagNatural,x,y,time,normalVector,varargin{:});
    ind = mesh.(edges)(4:3+nodesOnEdge,dof.(naturalName)(i));
    rhs.(domainVariable)(ind) = rhs.(domainVariable)(ind) + ...
        (edgeLength/2)*((wi.*fnxy)*phi')';
end
    
end  