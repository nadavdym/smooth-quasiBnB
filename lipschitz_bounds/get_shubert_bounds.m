function [c0,f,g,H,L1,L2,L3,gt_sol]=get_shubert_bounds()

%% get cube [-10,10]*[-10,10]
c0.h=[10;10];
c0.x=zeros(2,1);
gt_sol=nan;
filename='shubert_bound';
%% construct the function f

x=sym('x',[2,1]);
sum1 = 0; sum2 = 0;
for i = 1:5
 sum1 = sum1 + i.*cos((i+1).*x(1)+i);
 sum2 = sum2 + i.*cos((i+1).*x(2)+i);
end
f_symb = sum1.*sum2;


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