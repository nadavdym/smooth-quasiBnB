function [x,local_ub,qlb,fun_eval_num,lin_solve_num]=compute_lip_bounds(cube,f,L1)
%compute second order bounds

x=cube.x;
local_ub=f(x);
qlb=local_ub-L1*norm(cube.h);

fun_eval_num=1;
lin_solve_num=0;
end