function [u] = setStokesdarcySolution (w,dirichletValues,dof,domainStokes,domainDarcy,variableDarcy)

numberOfInternalStokes      = strcat(domainStokes,'_v1_ninternal');
numberOfPressureNodesStokes = strcat(domainStokes,'_nv2');

% Set Stokes solution
numberDofsStokes = dof.(numberOfInternalStokes) + dof.(numberOfPressureNodesStokes);
wStokes = w(1:numberDofsStokes,1);
% Values of Dirichlet boundary conditions
[u] = setStokesBC(wStokes,dirichletValues,dof,domainStokes);

% Set Darcy solution
wDarcy = w(numberDofsStokes+1:end,1);
% Values of Dirichlet boundary conditions
[u] = setLaplaceBC(wDarcy,dirichletValues,dof,domainDarcy,variableDarcy,u);

end