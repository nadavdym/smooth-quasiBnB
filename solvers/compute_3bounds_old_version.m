function [y,local_ub,qlb]=compute_3bounds(cube,f,grad,hess,L3)
%compute second order bounds
b=2.01;
a=8*L3/3;

x=cube.x;
diam=norm(cube.h);
g=grad(x);
H=hess(x);
[V,D]=eig(H);
Dvec=diag(D);
Dvec_proj=max(Dvec,zeros(size(Dvec)));
Dvec_new=Dvec_proj+a*diam/2;
y_temp=-diag(1./Dvec_new)*V'*g;
y=V*y_temp+x;
local_ub=f(y);
if norm(Dvec_proj-Dvec,'inf')>L3*diam || norm(y-x)>b*diam
    qlb=inf;
else
    delta=y-x;
    qlb=f(x)+g'*delta+0.5*delta'*V*diag(Dvec_new)*V'*delta-(2*L3/3+a/2)*diam^3;
end