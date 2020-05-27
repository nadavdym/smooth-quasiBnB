function [c0,f,g,H,L1,L2,L3,gt_sol]=get_gold_bounds()

%% get cube [-3,3]*[-2,2]
c0.h=[2;2];
c0.x=zeros(2,1);
gt_sol=[0 ; -1];
filename='gold_bound';
%% construct the function f

x=sym('x',[2,1]);
f_symb=(1+(x(1)+x(2)+1).^2.*(19-14.*x(1)+3.*x(1).^2-14.*x(2)+6.*x(1).*x(2)+3.*x(2).^2))...
.*(30+(2.*x(1)-3.*x(2)).^2.*(18-32.*x(1)+12.*x(1).^2+48.*x(2)-36.*x(1).*x(2)+27.*x(2).^2));
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