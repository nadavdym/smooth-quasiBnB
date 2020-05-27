function [minimizer,ub,qlb,fun_eval_num,lin_solve_num]=compute_alpha_bounds(cube,f,g,H,alpha_vec)
lb_cube=cube.x-cube.h;
ub_cube=cube.x+cube.h;
F=@(x) combine(x,f,g,alpha_vec,lb_cube,ub_cube);
display='off'; %'final' or 'off'
newH=@(y) H(y)+2*diag(alpha_vec);
[minimizer,qlb,fun_eval_num,lin_solve_num]=solve_on_cube(cube,F,newH,display);
ub=f(minimizer);



end

function [F,G] = combine(x,f,g,alpha_vec,lb_cube,ub_cube)
% Calculate objective f
F = f(x)+sum( alpha_vec.*(x-ub_cube).*(x-lb_cube) );

if nargout ==2 % gradient required
    G = g(x)+alpha_vec.*(x-lb_cube)+alpha_vec.*(x-ub_cube);
end
end

