function [mesh] = meshP1toP1B3 (mesh)

% meshP1toP1B3 computes the grid for P1-cubic bubble elements (MINI)
% starting from a P1 mesh.
%
% Code adapted from the original function P1toB1mesh2D.m by F. Saleri.

fem = int2str(mesh.t(4,1));
fem = strcat('d',fem,'_fem');
mesh.(fem) = 'B1';

mesh.t = [mesh.t(1:3,:); [1:size(mesh.t,2)]+size(mesh.p,2); mesh.t(4,:)];
for ie = 1:size(mesh.t,2)
    xb = sum(mesh.p(1,mesh.t(1:3,ie)))/3;
    yb = sum(mesh.p(2,mesh.t(1:3,ie)))/3;
    mesh.p = [mesh.p, [xb;yb]];
end

end
