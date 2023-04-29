function [u] = extractSolutionDarcyScalarStokesBjs (u,dof,domainStokes,domain)

dname = strcat('d',domainStokes);
internalS = strcat(dname,'_internal');
M = dof.(internalS)(1,2);
nInt = M - 2;

npnodes  = strcat(dname,'_v3');

if strcmp(domain,'1')     % Stokes domain
    u = u(1:nInt+dof.(npnodes));
elseif strcmp(domain,'2') % Darcy domain
    u = u(1+nInt+dof.(npnodes):end);
end

return