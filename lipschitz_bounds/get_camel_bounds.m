function [c0,f,g,H,L1,L2,L3,gt_sol]=get_camel_bounds()

%% get cube [-3,3]*[-2,2]
c0.h=[3;2];
c0.x=zeros(2,1);
x_star=[0.0898; -0.7126];
gt_sol=[x_star -x_star];
filename='camel_bound';
%% construct the function f


x=sym('x',[2,1]);
f_symb=(4-2.1.*x(1).^2+x(1).^4./3).*x(1).^2+x(1).*x(2)+(-4+4.*x(2).^2).*x(2).^2;
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
L=load(filename);
L=L.Expression1;
L1=L(1);
L2=L(2);
L3=L(3);









end