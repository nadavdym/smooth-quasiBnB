function [c0,f,g,H,L1,L2,L3,gt_sol]=get_branin_bounds()

%% get cube [-5,10]*[0,15]
c0.h=0.5*[15;15];
c0.x=0.5*[10-5;15];
gt_sol=[-pi, 12.275; pi 2.275; 9.42478, 2.475]';


M=[10;15]; %maximal magnitude of each coordinate

%% construct the function f

a=1;
b=5.1/(4*pi*pi);
c=5/pi;
d=6;
h=10;
ff=1/(8*pi);
x=sym('x',[2,1]);
f_symb=a.*(x(2)-b.*x(1).^2+c.*x(1)-d).^2+h.*(1-ff).*cos(x(1))+h;
%f_symb=x(1)^2+x(2)^2;
f_temp=matlabFunction(f_symb);
f=@(y) f_temp(y(1),y(2));
%% get gradient and hessian
g_symb=gradient(f_symb,x);
g_temp=matlabFunction(g_symb);
g=@(y) g_temp(y(1),y(2));
H_symb=hessian(f_symb,x);
H_temp=matlabFunction(H_symb);
H=@(y) H_temp(y(1),y(2));
%load bounds
L=load('branin_bound');
L=L.Expression1;
L1=L(1);
L2=L(2);
L3=L(3);









end