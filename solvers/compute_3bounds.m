function [y,local_ub,qlb,fun_eval_num,lin_solve_num]=compute_3bounds(cube,f,grad,hess,L3,tol)
%compute second order bounds
x0=cube.x;
fun_eval_num=1; %number of function evaluations
local_ub=f(x0); %this is default if newton iterations fail
y=x0;
d=length(x0);
diam=norm(cube.h);
g=grad(x0);
H=hess(x0);
e=eig(H);
lin_solve_num=1; %number of linear solves/eig decomposition
lambda_min=min(e);
lambda_max=max(e);
if (lambda_min<-L3*diam) 
    qlb=inf;
else
   lambda_plus=max(0,5*L3*diam-lambda_min);
   m=3*L3*diam; %minimial eigenvalue of new_f defined as f+lambda_plua/2*|x-x_0|^2 is at least m
   M=lambda_max+lambda_plus+L3*2*diam;
   x=x0;
   r_prev=diam;
   qlb=-inf;
   iter_count=0;
   while  (qlb==-inf)
      g=grad(x)+lambda_plus*(x-x0);
      H=hess(x)+lambda_plus*eye(d);
      x_prev=x;
      x=x-H\g;
      r_new=3*L3/(2*m)*r_prev^2;
      if (norm(x-x0)>r_new+diam) || norm(x-x_prev)>r_new+r_prev
          qlb=inf;
      else
          err=0.5*M*r_new^2; %possible error in minimizing new_f
          if err<tol
              qlb=f(x)-err-lambda_plus/2*diam^2;
              local_ub=f(x);
              y=x;
          end
      end
      iter_count=iter_count+1;
      lin_solve_num=lin_solve_num+1;
      r_prev=r_new;
   end
end

end



