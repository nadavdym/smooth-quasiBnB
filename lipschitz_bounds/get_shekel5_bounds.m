function [c0,f,g,H,L1,L2,L3,gt_sol]=get_shekel5_bounds()

%% get cube [-10,10]*[-10,10]
m=5;
c0.h=5*ones(4,1);
c0.x=5*ones(4,1);
gt_sol=4*ones(4,1);
filename='shekel5_bound';
%% construct the function f

x=sym('x',[4,1]);

a = [4, 1, 8, 6, 3, 2, 5, 8, 6, 7;
     4, 1, 8, 6, 7, 9, 5, 1, 2, 3.6;
     4, 1, 8, 6, 3, 2, 3, 8, 6, 7;
     4, 1, 8, 6, 7, 9, 3, 1, 2, 3.6];
c = [ 0.1, 0.2, 0.2, 0.4, 0.4, 0.6, 0.3, 0.7, 0.5, 0.5];
a=a(:,1:m);
c=c(1:m);
for i=1:m
 b = (x - a(:,i)).^2;
 d(i) = sum(b);
end
f_symb = -sum((c+d).^(-1));





f_temp=matlabFunction(f_symb);
f=@(y) f_temp(y(1),y(2),y(3),y(4));
%% get gradient and hessian
g_symb=gradient(f_symb,x);
g_temp=matlabFunction(g_symb);
g=@(y) g_temp(y(1),y(2),y(3),y(4));
H_symb=hessian(f_symb,x);
H_temp=matlabFunction(H_symb);
H=@(y) H_temp(y(1),y(2),y(3),y(4));
%load bounds
L=load(filename);
L=L.Expression1;
L1=L(1);
L2=L(2);
L3=L(3);









end