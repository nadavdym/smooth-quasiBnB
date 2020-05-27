



function [minimizer,fval,fun_eval_num,lin_solve_num]=solve_on_cube(cube,F,H,display)
%solve optimization problems on cube
if nargin==2
    display='off';
end
x0=cube.x;
lb=x0-cube.h;
ub=x0+cube.h;
A = [];
b = [];
Aeq = [];
beq = [];
nonlcon=[];
%options = optimoptions('fmincon','SpecifyObjectiveGradient',true,'Display',display);
hessianfcn=@(x,lambda) compute_hessian(x,lambda,H);
solver='sqp';
switch solver
    case 'interior-point'
        options = optimoptions('fmincon','Algorithm','interior-point',...
            'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
            'HessianFcn',hessianfcn,'Display',display,'SubproblemAlgorithm','factorization');
    case 'trust-region-reflective'
        options = optimoptions('fmincon','Algorithm','trust-region-reflective',...
            'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
            'HessianFcn',hessianfcn,'Display',display,'SubproblemAlgorithm','factorization');
    case 'sqp'
        options = optimoptions('fmincon','Algorithm','sqp',...
            'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
            'Display',display);
    case 'active-set'
        options = optimoptions('fmincon','Algorithm','active-set',...
            'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
            'Display',display);
    case 'sqp-legacy'
           options = optimoptions('fmincon','Algorithm','sqp-legacy',...
            'SpecifyObjectiveGradient',true,'SpecifyConstraintGradient',true,...
            'Display',display);
end
[minimizer,fval,exitflag,output] = fmincon(F,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

fun_eval_num=output.funcCount;
lin_solve_num=inf;

end

function Heval=compute_hessian(x,lambda,H)
Heval=H(x);
%Heval=eye(2);
%warning('bla');
end




