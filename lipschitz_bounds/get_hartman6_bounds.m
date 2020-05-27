function [c0,f,g,H,L1,L2,L3,gt_sol]=get_hartman6_bounds()

%% get cube [-10,10]*[-10,10]
c0.h=0.5*ones(6,1);
c0.x=0.5*ones(6,1);
gt_sol=[0.20169; 0.150011; 0.476874; 0.275332; 0.311652; 0.6573];
filename='hartman6_bound';
%% construct the function f

x=sym('x',[6,1]);

a = [10.00,  0.05,  3.00, 17.00;
     3.00, 10.00,  3.50,  8.00;
     17.00, 17.00,  1.70,  0.05;
     3.50,  0.10, 10.00, 10.00;
     1.70,  8.00, 17.00,  0.10;
     8.00, 14.00,  8.00, 14.00];
p = [0.1312, 0.2329, 0.2348, 0.4047;
     0.1696, 0.4135, 0.1451, 0.8828;
     0.5569, 0.8307, 0.3522, 0.8732;
     0.0124, 0.3736, 0.2883, 0.5743;
     0.8283, 0.1004, 0.3047, 0.1091;
     0.5886, 0.9991, 0.6650, 0.0381];
c = [1.0, 1.2, 3.0, 3.2];
for i=1:4
 d(i) = sum(a(:,i).*(x - p(:,i)).^2);
end

f_symb = -sum(c.*exp(-d)); 


f_temp=matlabFunction(f_symb);
f=@(y) f_temp(y(1),y(2),y(3),y(4),y(5),y(6));
%% get gradient and hessian
g_symb=gradient(f_symb,x);
g_temp=matlabFunction(g_symb);
g=@(y) g_temp(y(1),y(2),y(3),y(4),y(5),y(6));
H_symb=hessian(f_symb,x);
H_temp=matlabFunction(H_symb);
H=@(y) H_temp(y(1),y(2),y(3),y(4),y(5),y(6));
%load bounds
L=load(filename);
L=L.Expression1;
L1=L(1);
L2=L(2);
L3=L(3);









end