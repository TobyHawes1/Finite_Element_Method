function [u] = setLaplacePeriodicBC (u,dof,domain,variable)

domainVariable  = strcat(domain,'_',variable);
periodicBCnodes = strcat(domainVariable,'_periodicBC');
index1 = dof.(periodicBCnodes)(1,:);
index2 = dof.(periodicBCnodes)(2,:);
u.(domainVariable)(index2,1) = u.(domainVariable)(index1,1);

end