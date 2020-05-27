function [x,local_ub,qlb,fun_eval_num,lin_solve_num]=compute_2bounds_constrained(cube,c0,f,L2)

fun_eval_num=1;
lin_solve_num=0;
%compute second order bounds
if min(c0.h-cube.h)==0   %don't do anything 
    x=cube.x;
    local_ub=f(x);
    qlb=-inf;
else
    M=cube.x+cube.h;
    m=cube.x-cube.h;
    M0=c0.x+c0.h;
    m0=c0.x-c0.h;
    boolM=(M==M0);
    bool_m=(m0==m);
    x=boolM.*M+bool_m.*m+0.5*(1-bool_m).*(1-boolM).*(M+m);
    local_ub=f(x);
    coordinate_wise_max=cube.h+(boolM+bool_m).*cube.h;
    diam_squared=coordinate_wise_max'*coordinate_wise_max;
    qlb=local_ub-0.5*L2*diam_squared;
end

end