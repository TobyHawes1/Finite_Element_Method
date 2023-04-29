function [meshP2] = meshP1toP2 (meshP1)

% MeshP1toP2 modifies an exisiting mesh for P1 elements into a mesh for P2
% elements by adding the new degrees of freedom

meshP2.p = zeros(2,2*size(meshP1.p,2)+size(meshP1.t,2)-1);
meshP2.p(:,1:size(meshP1.p,2)) = meshP1.p;
meshP2.e = [meshP1.e(1:2,:); zeros(1,size(meshP1.e,2)); meshP1.e(3:end,:)];
meshP2.t = [meshP1.t(1:3,:); zeros(3,size(meshP1.t,2)); meshP1.t(4,:)];

P2_p_number = size(meshP1.p,2);

for i=1:size(meshP1.t,2)
    x_P1 = zeros(2,3);
    for k=1:3
        x_P1(:,k) = meshP1.p(:,meshP1.t(k,i));
    end
    x_P2 = [x_P1(:,1)+x_P1(:,2) x_P1(:,2)+x_P1(:,3) x_P1(:,1)+x_P1(:,3)]/2;
    % are these vertices already known?
    P2_number = zeros(1,3);
    for k=1:3
        xexists = find(meshP2.p(1,1:P2_p_number)==x_P2(1,k));
        yexists = find(meshP2.p(2,1:P2_p_number)==x_P2(2,k));
        exists = intersect(xexists,yexists);
        if isempty(exists)
            P2_p_number = P2_p_number + 1;
            P2_number(k) = P2_p_number;
            % vertices of the P2 mesh
            meshP2.p(:,P2_p_number) = x_P2(:,k);
        else
            P2_number(k) = exists;
        end
    end
    % insert the new nodes in the element structure
    meshP2.t(4:6,i) = P2_number';
end
%
for i=1:size(meshP2.e,2)
    x = mean(meshP2.p(1,meshP2.e(1:2,i)));
    y = mean(meshP2.p(2,meshP2.e(1:2,i)));
    testx = abs(meshP2.p(1,:)-x);
    testy = abs(meshP2.p(2,:)-y);
    index = zeros(1,length(testx));
    n = 0;
    for j=1:length(testx)
        if ( testx(j)<eps & testy(j)<eps )
            n = n + 1;
            index(n) = j;
        end
    end
    index = index(1,1:n);
    meshP2.e(3,i) = index;
end

if isfield(meshP1,'fem')
    meshP2.fem = meshP1.fem;
end
if isfield(meshP1,'d')
    meshP2.d = meshP1.d;
end

end