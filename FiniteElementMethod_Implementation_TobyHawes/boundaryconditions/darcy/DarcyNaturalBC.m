function [rhs] = DarcyNaturalBC (rhs,mesh,dof,domain,data,type,varargin)

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

dname        = strcat('d',domain);
fem          = strcat(dname,'_fem');
natural_name = strcat(dname,'_',type,'_3');
nodes        = strcat(dname,'_p');
edges        = strcat(dname,'_e');

% Load quadrature nodes and weights
if strcmp(mesh.(fem),'P11')
    [ti,wi,phi] = BasisOnQuad1D('P1',10);
    bind = 2;
    edgefem = [3;4];
elseif strcmp(mesh.(fem),'P21')
    [ti,wi,phi] = BasisOnQuad1D('P1',10);
    bind = 2;
    edgefem = [4;5];
elseif ( strcmp(mesh.(fem),'P22') | strcmp(mesh.(fem),'P12') )
    [ti,wi,phi] = BasisOnQuad1D('P2',10);
    bind = 3;
    edgefem = [4;5];
end

if strcmp(type,'neumann')
    dataflag = 22; % Neumann
elseif strcmp(type,'robin')
    dataflag = 23; % Robin
end

for i=1:length(dof.(natural_name))
    x = mesh.(nodes)(1,mesh.(edges)(1:2,dof.(natural_name)(i)));
    y = mesh.(nodes)(2,mesh.(edges)(1:2,dof.(natural_name)(i)));
    side_length = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
    nv = mesh.(edges)(edgefem,dof.(natural_name)(i));
    x = ((1-ti)*x(1)+(ti)*x(2));
    y = ((1-ti)*y(1)+(ti)*y(2));
    fnxy = feval(data,dataflag,x,y,nv,varargin{:});
    ind = mesh.(edges)(1:bind,dof.(natural_name)(i));
    rhs.(dname)(ind) = rhs.(dname)(ind) + side_length*((wi.*fnxy)*phi')';
end
    
return  