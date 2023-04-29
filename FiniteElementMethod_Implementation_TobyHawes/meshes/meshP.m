function [mesh] = meshP(p,e,t,k)

%k -> degree of polynomial used (1-linear, 2-quadratic, 3-cubic)
%To call linear case: [mesh] = meshP(p,e,t,1)
%To call quadratic case: [mesh] = meshP(p,e,t,2)
%To call cubic case (classical): [mesh] = meshP(p,e,t,3)
%To call cubic case (Gauss-Lobatto): [mesh] = meshP(p,e,t,31)

mesh = struct;

switch k
    case 1
        [mesh1] = meshP1(p,e,t);
        mesh = mesh1;
    case 2
        [mesh1] = meshP1(p,e,t);
        [mesh2] = meshP2(mesh1);
        mesh = mesh2;
        disp("Mesh2d-2")
    case 3
        [mesh1] = meshP1(p,e,t);
        %P1 -> number of nodes in mesh in P1 case. Needed for a calculation
        %in MeshP3.
        p1 = size(mesh1.p,2);
        [mesh2] = meshP2(mesh1);
        [mesh3] = meshP3(mesh2, p1);
        mesh = mesh3;
    case 31
        [mesh1] = meshP1(p,e,t);
        %P1 -> number of nodes in mesh in P1 case. Needed for a calculation
        %in MeshP3_2.
        p1 = size(mesh1.p,2);
        [mesh2] = meshP2(mesh1);
        [mesh3] = meshP3gl(mesh2, p1);
        mesh = mesh3;
end

