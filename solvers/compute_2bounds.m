function [x,local_ub,qlb,fun_eval_num,lin_solve_num]=compute_2bounds(cube,f,L2)
%compute second order upper and lower bound using the derivative

x=cube.x;
local_ub=f(x);
qlb=f(x)-0.5*L2*cube.h'*cube.h;
fun_eval_num=1;
lin_solve_num=0;


end