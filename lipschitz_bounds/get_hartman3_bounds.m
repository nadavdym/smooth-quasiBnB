function [c0,f,g,H,L1,L2,L3,gt_sol]=get_hartman3_bounds()

%% get cube [-10,10]*[-10,10]
c0.h=0.5*[1;1;1];
c0.x=0.5*ones(3,1);
gt_sol=[0.114614; 0.555649; 0.852547];
filename='hartman3_bound';
%% construct the function f

x=sym('x',[3,1]);
a = [3.0d0,  0.1d0,  3.0d0,  0.1d0;
     10.0d0, 10.0d0, 10.0d0, 10.0d0;
     30.0d0, 35.0d0, 30.0d0, 35.0d0];
p = [ 0.36890d0, 0.46990d0, 0.10910d0, 0.03815d0;
      0.11700d0, 0.43870d0, 0.87320d0, 0.57430d0;
      0.26730d0, 0.74700d0, 0.55470d0, 0.88280d0];
c = [1.0d0, 1.2d0, 3.0d0, 3.2d0];
for i=1:4
     d(i) = sum(a(:,i).*(x - p(:,i)).^2);
end
f_symb = -sum(c.*exp(-d)); 


f_temp=matlabFunction(f_symb);
f=@(y) f_temp(y(1),y(2),y(3));
%% get gradient and hessian
g_symb=gradient(f_symb,x);
g_temp=matlabFunction(g_symb);
g=@(y) g_temp(y(1),y(2),y(3));
H_symb=hessian(f_symb,x);
H_temp=matlabFunction(H_symb);
H=@(y) H_temp(y(1),y(2),y(3));
%load bounds
L=load(filename);
L=L.Expression1;
L1=L(1);
L2=L(2);
L3=L(3);









end