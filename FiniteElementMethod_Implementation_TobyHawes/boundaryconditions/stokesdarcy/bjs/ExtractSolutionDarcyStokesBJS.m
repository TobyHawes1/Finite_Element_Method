function [u] = ExtractSolutionDarcyStokesBJS (u,dof,domainS,domain)

dname = strcat('d',domainS);
nvnodes  = strcat(dname,'_ninternal_1');
npnodes  = strcat(dname,'_p');
vi2nodes = strcat(dname,'_internal_2');
nv2nodes = strcat(dname,'_ninternal_2');
if ~isfield(dof,vi2nodes)
    nv2nodes = nvnodes;
end

if strcmp(domain,'1')     % Stokes domain
    u = u(1:dof.(nvnodes)+dof.(nv2nodes)+dof.(npnodes));
elseif strcmp(domain,'2') % Darcy domain
    u = u(1+dof.(nvnodes)+dof.(nv2nodes)+dof.(npnodes):end);
end

return