function [mesh1] = meshP1(p,e,t)
% n -> number of columns in matrix e
n = size(e,2);

% Creating row vectors to include x and y values respectively of the 
% normal vector,
normx = zeros(1,n);
normy = zeros(1,n);

% Calculating normal vectors
for i = 1:n
    
    % Extract (x,y) coordinates of the points of each triangle stored in
    % matrix e
    x1y1 = p(:,e(1,i));
    x2y2 = p(:,e(2,i));
    
    % Calculate difference between x and y coordinates using the direction
    % of travel given in the matrix e
    dx = x2y2(1) - x1y1(1);
    dy = x2y2(2) - x1y1(2);
    
    % Normalize vector
    len_norm = sqrt( dx^2 + dy^2 );
    dx = dx / len_norm;
    dy = dy / len_norm;
    
    % Find out which direction the normal vectors should point, the sixth
    % row of matrix e indicates whether the left hand side (when travelling 
    % the direction of the vector (clockwise)) is outside the domain (0) or 
    % inside the domain (1)
    switch e(6,i)
        case 0
            normx(i) = -dy;
            normy(i) = dx;
        case 1
            normx(i) = dy;
            normy(i) = -dx;
    end
end 

%Rearrange and remove rows that are not needed from matrices e and t
e(3,:) = e(2,:);
e(2,:) = e(1,:);
e(1,:) = e(5,:);

e(4,:) = [];
e(4,:) = [];
e(4,:) = [];
e(4,:) = [];

t(4,:) = [];

% Add the x and y coordinates of normal vectors to matrix e
e = [e(1:1, :); normx; e(2:end, :)];
e = [e(1:2, :); normy; e(3:end, :)];

% Create the mesh and add the fields p, e, t
mesh1 = struct;

mesh1.p = p;
mesh1.e = e;
mesh1.t = t;
