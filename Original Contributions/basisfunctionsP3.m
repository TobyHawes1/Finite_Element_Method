% Basis functions expressed as:
% chi(xi,yi) = c1*xi^3 + c2*yi^3 + c3*xi^2*yi + c4*xi*yi^2 + ...
%              c5*xi^2 + c6*yi^2 + c7*xi*yi + c8*xi + c9*yi + c10
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)
% 4 (1/3,0)
% 5 (2/3,0)
% 6 (2/3,1/3)
% 7 (1/3,2/3)
% 8 (0,2/3)
% 9 (0,1/3)
% 6 (1/3,1/3)

xi = [0 1 0 1/3 2/3 2/3 1/3 0 0 1/3];
yi = [0 0 1 0 0 1/3 2/3 2/3 1/3 1/3];

A = zeros(10,10);

for i = 1:10
    A(i,:) = [xi(i)^3, yi(i)^3, xi(i)^2*yi(i), xi(i)*yi(i)^2, xi(i)^2, yi(i)^2, xi(i)*yi(i), xi(i), yi(i), 1];
end

b = eye(10,10);

% Coefficients of the basis functions
c = A\b;

P1 = [xi(1,1),yi(1,1)];
P2 = [xi(1,2),yi(1,2)];
P3 = [xi(1,3),yi(1,3)];
% Create triangular mesh of triangla P1,P2,P3
n=20;
[I,J]=meshgrid(linspace(0,1,n));
keep=I+J <= 1+eps;
IJ = [I(keep),J(keep)];
K = 1-sum(IJ,2);
F=delaunay(IJ);
IJK = [IJ,K];
XY=IJK*[P1;P2;P3];
x=XY(:,1);
y=XY(:,2);
close all
for j=1:size(A,1)
    phi = @(x,y) c(1,j)*x.^3 + c(2,j)*y.^3 + c(3,j)*x.^2.*y + c(4,j)*x.*y.^2 + c(5,j)*x.^2 + c(6,j)*y.^2 + c(7,j)*x.*y + c(8,j)*x + c(9,j)*y + c(10,j);
    figure(j)
    trisurf(F,x,y,phi(x,y),'FaceColor',(192/255)*[1 1 1],'EdgeColor',(169/255)*[1 1 1])
end

syms x y
chi = sym(zeros(10,1));
dchidx = sym(zeros(10,1));
dchidy = sym(zeros(10,1));
for j=1:10
    chi(j,1) = sym(c(1,j),'r')*x^3 + ...
               sym(c(2,j),'r')*y^3 + ...
               sym(c(3,j),'r')*x^2*y + ...
               sym(c(4,j),'r')*x*y^2 + ...
               sym(c(5,j),'r')*x^2 + ...
               sym(c(6,j),'r')*y^2 + ...
               sym(c(7,j),'r')*x*y + ...
               sym(c(8,j),'r')*x + ...
               sym(c(9,j),'r')*y + ...
               sym(c(10,j),'r');
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end

symDERXX = sym(zeros(10,10));
symDERYY = sym(zeros(10,10));
symDERXY = sym(zeros(10,10));
symMASS = sym(zeros(10,10));
for i = 1:10
    for j = 1:10
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),y,0,-x+1),x,0,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),y,0,-x+1),x,0,1);
    end
end