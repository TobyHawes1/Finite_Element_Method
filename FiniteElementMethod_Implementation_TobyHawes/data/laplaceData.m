function [varargout] = laplaceData (whichData,varargin)

% WhichData is a flag with the following correspondences:
% 10 - diffusion coefficient
% 11 - force
% 20 - type of boundary
% 21 - Dirichlet boundary condition
% 22 - Neumann boundary condition
% 23 - Robin boundary condition
% 24 - Robin coefficient
% 30 - exact solution

% - diffusion * Laplacian(u) = force

switch whichData
    case 10
        % DIFFUSION coefficient of the Laplace equation
        x = varargin{1};
        y = varargin{2};
        diffusion = 1 + 0*(x+y);
        varargout(1) = {diffusion};
        %
    case 11
        % FORCE (scalar function!) of the Laplace equation
        x = varargin{1};
        y = varargin{2};
        force = 8*pi^2*sin(2*pi*x).*cos(2*pi*y);
        %force = 0*x;
        varargout(1) = {force};
        %
    case 20
        % TYPE of BOUNDARY
        boundary.d1_DirichletEdges = [1 2 3 4];%[1 2 3 4 5];
        boundary.d1_NeumannEdges   = [];
        boundary.d1_RobinEdges     = [];
        varargout(1) = {boundary};
        %
    case 21
        % DIRICHLET boundary conditions
        x = varargin{1};
        y = varargin{2};
        dirichlet = laplaceData(30,x,y);
        varargout(1) = {dirichlet};
        %
    case 22
        % NEUMANN boundary conditions
        x = varargin{1};
        y = varargin{2};
        normalvector = varargin{4};
        [u,ux,uy] = laplaceData(30,x,y);
        neumann = laplaceData(10,x,y).*(ux*normalvector(1)+uy*normalvector(2));
        varargout(1) = {neumann};
        %
    case 23
        % ROBIN boundary conditions
        x = varargin{1};
        y = varargin{2};
        normalvector = varargin{4};
        [u,ux,uy] = laplaceData(30,x,y);
        robin = laplaceData(10,x,y).*(ux*normalvector(1)+uy*normalvector(2)) + laplaceData(24)*u;
        varargout(1) = {robin};
        %
    case 24
        % ROBIN coefficient
        robincoefficient = 1;
        varargout(1) = {robincoefficient};
        %
    case 30
        % EXACT SOLUTION
        x = varargin{1};
        y = varargin{2};
        u = sin(2*pi*x).*cos(2*pi*y);
        varargout(1) = {u};
        if nargout>1
            dudx =  2*pi*cos(2*pi*x).*cos(2*pi*y);
            dudy = -2*pi*sin(2*pi*x).*sin(2*pi*y);
            varargout(2) = {dudx};
            varargout(3) = {dudy};
        end
        %
end

end