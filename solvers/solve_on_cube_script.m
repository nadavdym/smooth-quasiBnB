function solve_on_cube_script


x0 = [1.1,1.1];
cube.x=x0;
cube.h=[0.3 0.3];
a=1; %this doens't do anything
F=@(x) rosenbrockwithgrad(x,a);
display='off'; %'final' or 'off'
[minimizer,fval]=solve_on_cube(cube,F,display);
end
% lb = [0.8, 0.8];
% ub = [1.2, 1.2];
% 
% A = [];
% b = [];
% Aeq = [];
% beq = [];
% nonlcon=[];
% fun=@rosenbrockwithgrad;
% options = optimoptions('fmincon','SpecifyObjectiveGradient',true);
% x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

% function solve_on_cube_script
% d=2;
% a=zeros(d,1);
% theta=2*pi;
% delta=10^-1;
% x_0=ones(d,1);
% x_0(1)=0;
% f=@(x) a'*(1-cos(theta*x))+delta*(x-x_0)'*(x-x_0);
% c0.x=rand(d,1);
% c0.h=2*ones(d,1);
% g=@(x) theta*a.*sin(theta*x)+2*delta*(x-x_0);
% H=@(x) compute_Hessian(x,a,theta,delta);
% [x,val]=solve_on_cube(c0,f,g,H);
% end
% 
% function H=compute_Hessian(x,a,theta,delta)
% Diag=theta^2*a.*cos(theta*x)+2*delta;
% H=diag(Diag);
% end
% 
% function [f,g]=join(f,g)
% f=f;
% if nargout>1
% g=g;
% end
% end

function [f,g] = rosenbrockwithgrad(x,a)
% Calculate objective f
f = 100*(x(2) - x(1)^2)^2 + (1-x(1))^2+a;

if nargout > 1 % gradient required
    g = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1));
        200*(x(2)-x(1)^2)];
end
end




