function [rhs] = DarcyStokesBJSCouplingBC (A,matdata,rhs,dof,dirichlet,domains)

stokesdomain = strcat('d',domains{1});
stkdomain    = str2double(domains{1});
stokesvalues = strcat(stokesdomain,'_2');
darcydomain  = strcat('d',domains{2});
drcdomain    = str2double(domains{2});
darcyvalues  = strcat(darcydomain,'_3');

% Correction to the Stokes problem due to the coupling
unodesS = strcat(stokesdomain,'_u');
unodesD = strcat(darcydomain,'_p');
gnodes  = strcat(stokesdomain,'_gamma_2');
if ~isfield(dof,gnodes)
    gnodes = strcat(stokesdomain,'_gamma_1');
end
dnodes = strcat(darcydomain,'_dirichlet_3');
temp = sparse(dof.(unodesD),1);
temp(dof.(dnodes),1) = dirichlet.(darcyvalues);
%dof.(dnodes) = intersect(dof.gamma,dof.(dnodes));
temp = temp(dof.(dnodes),1);
%
if ~isempty(temp)
    rhs.(stokesdomain)(dof.(gnodes)(stkdomain,:)+dof.(unodesS),1) = ...
        rhs.(stokesdomain)(dof.(gnodes)(stkdomain,:)+dof.(unodesS),1) + ...
        A.(matdata.MGamma)(dof.(gnodes)(stkdomain,:),dof.(dnodes))*temp;
end

% Correction to the Darcy problem due to the coupling
gnodes = strcat(darcydomain,'_gamma_3');
dnodes = strcat(stokesdomain,'_dirichlet_2');
if ~isfield(dof,dnodes)
    dnodes = strcat(stokesdomain,'_dirichlet_1');
end
temp = sparse(dof.(unodesS),1);
temp(dof.(dnodes),1) = dirichlet.(stokesvalues);
%dof.(dnodes) = intersect(dof.gamma,dof.(dnodes));
temp = temp(dof.(dnodes),1);
%
rhs.(darcydomain)(dof.(gnodes)(drcdomain,:),1) = ...
    rhs.(darcydomain)(dof.(gnodes)(drcdomain,:),1) - ...
    ((A.(matdata.MGamma)(dof.(dnodes),dof.(gnodes)(drcdomain,:)))')*temp;

return