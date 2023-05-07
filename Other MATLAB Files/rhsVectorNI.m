function [rhs] = rhsVectorNI (rhs,rhsname,rhstype,mesh,dof,domainVariable,fem,...
    numberOfNodes,elements,vind,data,dataflag,dcdx,dcdy,dedx,dedy,aree,time,varargin)

degree = 10;
switch fem
    case {'Q1','Q2','Q3','Q3gl'}
        degree = 5;
end

[csi,eta,w,phi,dphix,dphiy] = basisOnQuad2D(fem,degree);

% Assemble the right-hand side WITH numerical integration
switch rhstype
    case {1,3,4,9}
        % Scalar case
        rhs.(rhsname) = zeros(dof.(numberOfNodes),1);
        for ie = 1:size(mesh.(elements),2)
            ind = (mesh.(elements)(1:vind,ie))';
            rhsloc = localRhsVectorNI(data,dataflag,mesh,aree(ie),domainVariable,ind,...
                fem,degree,dcdx(ie),dcdy(ie),dedx(ie),dedy(ie),rhstype,time,...
                csi,eta,w,phi,dphix,dphiy);
            rhs.(rhsname)(ind) = rhs.(rhsname)(ind) + rhsloc(1:vind,1);
        end
    case {2,5,6,7,8}
        % Vector case
        rhs.(rhsname) = zeros(2*dof.(numberOfNodes),1);
        for ie = 1:size(mesh.(elements),2)
            ind1 = (mesh.(elements)(1:vind,ie))';
            ind2 = ind1 + dof.(numberOfNodes);
            rhsloc = localRhsVectorNI(data,dataflag,mesh,aree(ie),domainVariable,ind1,...
                fem,degree,dcdx(ie),dcdy(ie),dedx(ie),dedy(ie),rhstype,time,...
                csi,eta,w,phi,dphix,dphiy);
            rhs.(rhsname)(ind1) = rhs.(rhsname)(ind1) + rhsloc(1:vind,1);
            rhs.(rhsname)(ind2) = rhs.(rhsname)(ind2) + rhsloc(vind+1:2*vind,1);
        end
end
        
end