function [rhs,dirichlet,dof] = DarcyMixedDirichletBC (A,matdata,rhs,...
    mesh,dof,domain,data,settings,dirichlet,varargin)

% DarcyMixedDirichletBC imposes Dirichlet boundary conditions for the Darcy
% problem in mixed form.

dname = strcat('d',domain);
nameu = strcat(dname,'_1_2');
namep = strcat(dname,'_3');
internal        = strcat(dname,'_internal');
dirichlet_name1 = strcat(dname,'_dirichlet_1');
dirichlet_name2 = strcat(dname,'_dirichlet_2');

unodes = strcat(dname,'_u');
nodes  = strcat(dname,'_p');
diri1  = strcat(dname,'_1');
diri2  = strcat(dname,'_2');

N2 = dof.(internal)(1,2);

if settings.darcymixed<4
    % Impose in essential way the boundary condition on the normal velocity
    functionflag = 21;
    % Data "21" corresponds to Dirichlet bc
    % Dirichlet bc on the first component of the velocity
    if ~isempty(dof.(dirichlet_name1))
        %
        component = 1;
        x = mesh.(nodes)(1,abs(dof.(dirichlet_name1)));
        y = mesh.(nodes)(2,abs(dof.(dirichlet_name1)));
        dirichlet.(diri1) = feval(data,functionflag,x,y,component,varargin{:})';
        dirichlet.(diri1) = sign(dof.(dirichlet_name1)').*dirichlet.(diri1);
        dof.(dirichlet_name1) = abs(dof.(dirichlet_name1));
        % Momentum equation
        rhs.(nameu)(dof.(internal)(3:N2)) = ...
            rhs.(nameu)(dof.(internal)(3:N2)) - ...
            A.(matdata.Mass11)(dof.(internal)(3:N2),dof.(dirichlet_name1))*dirichlet.(diri1);
        % Continuity equation
        rhs.(namep) = rhs.(namep) - ...
            (A.(matdata.G)(dof.(dirichlet_name1),:)')*dirichlet.(diri1);
    end
    % Dirichlet bc on the second component of the velocity
    if ~isempty(dof.(dirichlet_name2))
        %
        component = 2;
        x = mesh.(nodes)(1,abs(dof.(dirichlet_name2)));
        y = mesh.(nodes)(2,abs(dof.(dirichlet_name2)));
        dirichlet.(diri2) = feval(data,functionflag,x,y,component,varargin{:})';
        dirichlet.(diri2) = sign(dof.(dirichlet_name2)').*dirichlet.(diri2);
        dof.(dirichlet_name2) = abs(dof.(dirichlet_name2));
        % Momentum equation
        rhs.(nameu)(dof.(internal)(3:N2)) = ...
            rhs.(nameu)(dof.(internal)(3:N2)) - ...
            A.(matdata.Mass11)(dof.(internal)(3:N2),...
            dof.(dirichlet_name2)+dof.(unodes))*dirichlet.(diri2);
        % Continuity equation
        rhs.(namep) = rhs.(namep) - ...
            (A.(matdata.G)(dof.(dirichlet_name2)+dof.(unodes),:)')*dirichlet.(diri2);
    end
elseif settings.darcymixed>3
    % Impose in essential way the boundary condition on the pressure
    functionflag = 22;
    dirichlet_name1 = strcat(dname,'_dirichlet_3');
    diri1           = strcat(dname,'_3');
    internal        = strcat(dname,'_internal');
    N1 = dof.(internal)(1,1);
    x = mesh.(nodes)(1,dof.(dirichlet_name1));
    y = mesh.(nodes)(2,dof.(dirichlet_name1));
    dirichlet.(diri1) = feval(data,functionflag,x,y,varargin{:})';
    % Momentum equation
    rhs.(nameu) = rhs.(nameu) - ...
        A.(matdata.G)(:,dof.(dirichlet_name1))*dirichlet.(diri1);
    % Continuity equation
    rhs.(namep)(dof.(internal)(1,2:N1),1) = rhs.(namep)(dof.(internal)(1,2:N1),1) - ...
        A.(matdata.StabQP)(dof.(internal)(1,2:N1),dof.(dirichlet_name1))*dirichlet.(diri1);
end

return