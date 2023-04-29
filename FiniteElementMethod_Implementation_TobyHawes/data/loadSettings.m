function [settings] = loadSettings ()

% settings.solver : solution methods
%                   'Direct' - solution by factorization (Cholesky)
%                   'PCG'    - CG with preconditioner from an incomplete
%                              factorization of the stiffness matrix
%                   'CG'     - CG without preconditioner
%
% settings.NI     : numerical integration
%                   1 - assembles the stiffness matrix using numerical
%                       integration
%                   0 - assembles the stiffness matrix without numerical
%                       integration
%
% setting.plot    : 0 - does not plot the computed solution
%                   1 - plots the computed solution
%
% settings.error  : 0 - does not compute any error
%                   1 - computes the error with respect to a known solution

settings.solver = 'Direct';

settings.NI     = 0;

settings.plot      = 1;
settings.plotExact = 1;

settings.error  = 1;

end