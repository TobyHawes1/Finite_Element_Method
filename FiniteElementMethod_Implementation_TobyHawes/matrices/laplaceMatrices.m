function [A,rhs,matdata,varargout] = laplaceMatrices (mesh,dof,domain,variable,...
    time,data,matdata,whichMatrix,varargin)

% LaplaceMatrices assembles the matrices and the rhs required to solve
% a scalar Laplace problem.

% In case A and/or rhs are already defined:
if ~isempty(varargin)
    A = varargin{1};
    if length(varargin)>1
        rhs = varargin{2};
    else
        rhs = {};
    end
else
    A   = {};
    rhs = {};
end

domainVariable = strcat(domain,'_',variable);
elements = strcat(domainVariable,'_t');
nodes    = strcat(domainVariable,'_p');
fem      = strcat(domainVariable,'_fem');
numberOfNodes = strcat(domain,'_n',variable);
fem = mesh.(fem);
%fem = fem(end-1:end);

switch fem
    case {'P1'}
        % Linear Lagrangian elements on triangles
        femflag = 0;
        factor  = 0.5;
        vind    = 3;
    case {'P2'}
        % Quadratic Lagrangian elements on triangles
        femflag = 0;
        factor  = 0.5;
        vind    = 6;
    case {'P3','P3gl'}
        % Cubic Lagrangian elements on triangles
        femflag = 0;
        factor  = 0.5;
        vind    = 10;
    case {'Q1'}
        % Linear Lagrangian elements on rectangles
        femflag = 1;
        factor  = 1;
        vind    = 4;
    case {'Q2'}
        % Quadratic Lagrangian elements on rectangles
        femflag = 1;
        factor  = 1;
        vind    = 9;
    case {'Q3','Q3gl'}
        % Cubic Lagrangian elements on rectangles
        femflag = 1;
        factor  = 1;
        vind    = 16;
end

[aree,dcdx,dcdy,dedx,dedy] = geoTrasf2D(mesh.(nodes),mesh.(elements),femflag);
varargout(1) = {factor*aree};

% -------------------------------------------
% STIFFNESS matrix  (K grad u, grad v)_\Omega
% -------------------------------------------
if whichMatrix(1)~=0
    matdata.C = strcat(domain,'_C');
    dataflag = 10;
    if whichMatrix(1)>0
        % Assemble stiffness matrix WITHOUT numerical integration
        [A] = gradgradMatrix(A,matdata.C,whichMatrix(1),mesh,dof,domain,fem,...
            numberOfNodes,nodes,elements,vind,data,dataflag,dcdx,dcdy,dedx,dedy,aree);
    elseif whichMatrix(1)<0
        % Assemble stiffness matrix WITH numerical integration
        [A] = gradgradMatrixNI(A,matdata.C,mesh,dof,domainVariable,fem,...
        numberOfNodes,elements,vind,data,dataflag,dcdx,dcdy,dedx,dedy,aree);
    end
end

% -------------------------------------------
% MASS matrix  (u, v)_\Omega
% -------------------------------------------
if whichMatrix(2)~=0
    matdata.M = strcat(domain,'_M');
    %dataflag = 10;
    %A.(matdata.C) = zeros(dof.(unodes),dof.(unodes));
    if whichMatrix(2)>0
        % Assemble stiffness matrix WITHOUT numerical integration
        [A] = massMatrix (A,matdata.M,whichMatrix(2),mesh,dof,fem,...
            numberOfNodes,nodes,elements,vind,data,[],aree);
        %
    elseif whichMatrix(2)<0
        % Assemble mass matrix WITH numerical integration
        % Not available yet.
    end
end

% ----------------
% RHS (f,v)_\Omega
% ----------------
if whichMatrix(3)>0
    dataflag = 11;
    % Assemble the right-hand side WITH numerical integration
    [rhs] = rhsVectorNI(rhs,domainVariable,1,mesh,dof,domainVariable,fem,...
        numberOfNodes,elements,vind,data,dataflag,dcdx,dcdy,dedx,dedy,aree,time);
end

% ------------------
% RHS (uD',v)_\Omega
% ------------------
if whichMatrix(4)>0
    dataflag = 41;
    rhsName = strcat(domain,'_DtDirichlet');
    % Assemble the right-hand side WITH numerical integration
    [rhs] = rhsVectorNI(rhs,rhsName,1,mesh,dof,domain,fem,nodes,...
        numberOfNodes,elements,vind,data,dataflag,dcdx,dcdy,dedx,dedy,aree,time);
end

end