% Basis functions expressed as:
% chi(xi,yi) = c1*xi + c2*yi + c3
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)

xi = [0 1 0];
yi = [0 0 1];

A = zeros(3,3);

for i = 1:3
    A(i,:) = [xi(i), yi(i), 1];
end

b = eye(3,3);

% Coefficients of the basis functions
c = A\b;

% plot the functions
close all
P1 = [xi(1,1),yi(1,1)];
P2 = [xi(1,2),yi(1,2)];
P3 = [xi(1,3),yi(1,3)];

% Create triangular mesh of triangle P1,P2,P3
n=10;
[I,J]=meshgrid(linspace(0,1,n));
keep=I+J <= 1+eps;
IJ = [I(keep),J(keep)];
K = 1-sum(IJ,2);
F=delaunay(IJ);
IJK = [IJ,K];
XY=IJK*[P1;P2;P3];
x=XY(:,1);
y=XY(:,2);

for j=1:size(A,1)
    phi = @(x,y) c(1,j)*x + c(2,j)*y + c(3,j);
    figure(j)
    
    trisurf(F,x,y,phi(x,y),'FaceColor',(192/255)*[1 1 1],'EdgeColor',(169/255)*[1 1 1])
    hold on
    plot(0,0,'r*')
    plot(0,1,'r*')
    plot(1,0,'r*')
end

syms x y
chi = sym(zeros(3,1));
dchidx = sym(zeros(3,1));
dchidy = sym(zeros(3,1));
for j=1:3
    chi(j,1) = sym(c(1,j),'r')*x + ...
               sym(c(2,j),'r')*y + ...
               sym(c(3,j),'r');
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end

symDERXX = sym(zeros(3,3));
symDERYY = sym(zeros(3,3));
symDERXY = sym(zeros(3,3));
symMASS = sym(zeros(3,3));
for i = 1:3
    for j = 1:3
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),y,0,-x+1),x,0,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),y,0,-x+1),x,0,1);
    end
end
