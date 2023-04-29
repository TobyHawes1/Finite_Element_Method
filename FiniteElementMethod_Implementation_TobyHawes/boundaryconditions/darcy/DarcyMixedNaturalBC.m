function [rhs] = DarcyMixedNaturalBC ...
    (rhs,mesh,dof,domain,data,type,settings,varargin)

% DarcyMixedNaturalBC imposes natural boundary conditions to the Darcy
%  problem in mixed form.
%
% 'type' is a string which can take the values 'neumann' or 'robin'.
%
% In case of NEUMANN conditions, we need to define the pressure p on the
% boundary.
%

dname  = strcat('d',domain);
fem    = strcat(dname,'_fem');
unodes = strcat(dname,'_u');
nodes  = strcat(dname,'_p');
edges  = strcat(dname,'_e');
nameu  = strcat(dname,'_1_2');
namep  = strcat(dname,'_3');

if strcmp(type,'neumann')
    functionflag = 22; % Neumann
elseif strcmp(type,'robin')
    functionflag = 23; % Robin
end

% Load quadrature nodes and weights
switch mesh.(fem)
    case 'P11'
        vfem    = 'P1';
        pfem    = 'P1';
        vbind   = 2;
        pbind   = 2;
        edgefem = [3;4];
    case 'P12'
        vfem    = 'P1';
        pfem    = 'P2';
        vbind   = 2;
        pbind   = 3;
        edgefem = [4;5];
    case 'P21'
        vfem    = 'P2';
        pfem    = 'P1';
        vbind   = 3;
        pbind   = 2;
        edgefem = [4;5];
    case 'P22'
        vfem    = 'P2';
        pfem    = 'P2';
        vbind   = 3;
        pbind   = 3;
        edgefem = [4;5];
end

if settings.darcymixed<4
    [ti,wi,phi] = BasisOnQuad1D(vfem,10);
    % Impose in natural way the boundary condition on the pressure
    natural_name = strcat(dname,'_',type,'_1');
    for ie=1:length(dof.(natural_name))
        x = mesh.(nodes)(1,mesh.(edges)(1:2,dof.(natural_name)(ie)));
        y = mesh.(nodes)(2,mesh.(edges)(1:2,dof.(natural_name)(ie)));
        side_length = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
        nv = mesh.(edges)(edgefem,dof.(natural_name)(ie));
        x = ((1-ti)*x(1)+(ti)*x(2));
        y = ((1-ti)*y(1)+(ti)*y(2));
        fnxy = feval(data,functionflag,x,y,[],nv,varargin{:});
        fnxy = (nv(1)+nv(2))*fnxy;
        ind = mesh.(edges)(1:vbind,dof.(natural_name)(ie));
        if nv(2)~=0
            ind = ind + dof.(unodes);
        end
        rhs.(nameu)(ind) = rhs.(nameu)(ind) - side_length*((wi.*fnxy)*phi')';
    end
elseif settings.darcymixed>3
    [ti,wi,phi] = BasisOnQuad1D(pfem,10);
    % Impose in natural way the boundary condition on the normal velocity
    functionflag = 21;
    natural_name = strcat(dname,'_',type,'_1');
    for ie=1:length(dof.(natural_name))
        x = mesh.(nodes)(1,mesh.(edges)(1:2,dof.(natural_name)(ie)));
        y = mesh.(nodes)(2,mesh.(edges)(1:2,dof.(natural_name)(ie)));
        side_length = sqrt((x(2)-x(1))^2+(y(2)-y(1))^2);
        x = ((1-ti)*x(1)+(ti)*x(2));
        y = ((1-ti)*y(1)+(ti)*y(2));
        fnxy = feval(data,functionflag,x,y,[],varargin{:});
        ind = mesh.(edges)(1:pbind,dof.(natural_name)(ie));
        rhs.(namep)(ind) = rhs.(namep)(ind) + side_length*((wi.*fnxy)*phi')';
    end
end

return