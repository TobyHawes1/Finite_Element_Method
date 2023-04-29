function [rhs] = darcyScalarNaturalBC (rhs,mesh,dof,domain,data,type,varargin)

% DarcyNaturalBC imposes natural boundary conditions to the Poisson
%    equations. The boundary conditions can be of Neumann or Robin type.
%
% 'type' is a string which can take the values 'neumann' or 'robin'.
%
% In case of NEUMANN conditions, we need to define the conormal derivative
% \nu \nabla u \cdot \mathbf{n}
%
% In case of ROBIN conditions, we need to define the function
% \nu \nabla u \cdot \mathbf{n} + a * u

dName       = strcat('d',domain);
fem         = strcat(dName,'_fem');
naturalName = strcat(dName,'_',type,'_v3');
nodes       = strcat(dName,'_p');
edges       = strcat(dName,'_e');
name3       = strcat(dName,'_v3');

% Load quadrature nodes and weights
switch mesh.(fem)
    case 'P1P1'
        [ti,wi,phi] = basisOnQuad1D('P1',10);
        bind = 2;
        edgefem = [3;4];
    case 'P2P1'
        [ti,wi,phi] = basisOnQuad1D('P1',10);
        bind = 2;
        edgefem = [4;5];
    case {'P2P2','P1P2'}
        [ti,wi,phi] = basisOnQuad1D('P2',10);
        bind = 3;
        edgefem = [4;5];
    case 'Q1Q1'
        [ti,wi,phi] = basisOnQuad1D('P1',10);
        bind = 2;
        edgefem = [3;4];
    case 'Q1Q2'
        [ti,wi,phi] = basisOnQuad1D('P2',10);
        bind = 3;
        edgefem = [4;5];
    case 'Q2Q1'
        [ti,wi,phi] = basisOnQuad1D('P1',10);
        bind = 2;
        edgefem = [4;5];
    case 'Q2Q2'
        [ti,wi,phi] = basisOnQuad1D('P2',10);
        bind = 3;
        edgefem = [4;5];
end

if strcmp(type,'neumann')
    dataFlag = 21; % Neumann
elseif strcmp(type,'robin')
    dataFlag = 23; % Robin
end

for i=1:length(dof.(naturalName))
    x = mesh.(nodes)(1,mesh.(edges)(1:2,dof.(naturalName)(i)));
    y = mesh.(nodes)(2,mesh.(edges)(1:2,dof.(naturalName)(i)));
    side_length = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    nv = mesh.(edges)(edgefem,dof.(naturalName)(i));
    x = ((1-ti)*x(1)+(ti)*x(2));
    y = ((1-ti)*y(1)+(ti)*y(2));
    fnxy = feval(data,dataFlag,x,y,nv,varargin{:});
    ind = mesh.(edges)(1:bind,dof.(naturalName)(i));
    rhs.(name3)(ind) = rhs.(name3)(ind) - side_length*((wi.*fnxy)*phi')';
end
    
end