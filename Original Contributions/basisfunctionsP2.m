% Basis functions expressed as:
% chi(xi,yi) = c1*xi^2 + c2*yi^2 + c3*xi*yi + c4*xi + c5*yi + c6
% Numbering of the nodes:
% 1 (0,0)
% 2 (1,0)
% 3 (0,1)
% 4 (1/2,0)
% 5 (1/2,1/2)
% 6 (0,1/2)

xi = [0 1 0 1/2 1/2 0];
yi = [0 0 1 0 1/2 1/2];

A = zeros(6,6);

for i = 1:6
    A(i,:) = [xi(i)^2, yi(i)^2, xi(i)*yi(i), xi(i), yi(i), 1];
end

b = eye(6,6);

% Coefficients of the basis functions
c = A\b;


P1 = [xi(1,1),yi(1,1)];
P2 = [xi(1,2),yi(1,2)];
P3 = [xi(1,3),yi(1,3)];
% Create triangular mesh of triangle P1,P2,P3
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
    phi = @(x,y) c(1,j)*x.^2 + c(2,j)*y.^2 + c(3,j)*x.*y + c(4,j)*x + c(5,j)*y + c(6,j);
    figure(j)
    trisurf(F,x,y,phi(x,y),'FaceColor',(192/255)*[1 1 1],'EdgeColor',(169/255)*[1 1 1])
    hold on
    plot(0,0,'r*')
    plot(0,1,'r*')
    plot(1,0,'r*')
    plot(1/2,0,'r*')
    plot(1/2,1/2,'r*')
    plot(0,1/2,'r*')
end

syms x y
chi = sym(zeros(6,1));
dchidx = sym(zeros(6,1));
dchidy = sym(zeros(6,1));
for j=1:6
    chi(j,1) = sym(c(1,j),'r')*x*x + ...
               sym(c(2,j),'r')*y*y + ...
               sym(c(3,j),'r')*x*y + ...
               sym(c(4,j),'r')*x + ...
               sym(c(5,j),'r')*y + ...
               sym(c(6,j),'r');
    dchidx(j,1) = diff(chi(j,1),x);
    dchidy(j,1) = diff(chi(j,1),y);
end

symDERXX = sym(zeros(6,6));
symDERYY = sym(zeros(6,6));
symDERXY = sym(zeros(6,6));
symMASS = sym(zeros(6,6));
for i = 1:6
    for j = 1:6
        symDERXX(i,j) = int(int(dchidx(i,1)*dchidx(j,1),y,0,-x+1),x,0,1);
        symDERYY(i,j) = int(int(dchidy(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symDERXY(i,j) = int(int(dchidx(i,1)*dchidy(j,1),y,0,-x+1),x,0,1);
        symMASS(i,j)  = int(int(chi(i,1)*chi(j,1),y,0,-x+1),x,0,1);
    end
end