function [c0,f,g,H,L1,L2,L3,gt]=get_cosine_exp_bounds(d,alpha,delta,x0)
%d=2
%alpha=10*ones(d,1);
theta=2*pi;


gt=zeros(d,1);

f=@(x) alpha'*(1-cos(theta*x))+delta*x'*x;

if nargin>3
    c0.x=x0;
else
c0.x=zeros(d,1);
end
c0.h=5.12*ones(d,1);

L1=theta*norm(alpha)+2*abs(delta)*norm(c0.h);
L2=theta^2*norm(alpha)+2*abs(delta);
L3=norm(alpha)*theta^3;

g=@(x) theta*alpha.*sin(theta*x)+2*delta*(x);
%H=@(x) compute_Hessian(x,alpha,theta,delta);

H=@(x) diag(theta^2*alpha.*cos(theta*x))+2*delta*eye(2);
%plot function
doDisplay=true;
if doDisplay && d==2
    [x,y]=meshgrid(-2:0.1:2);
    z=inf(size(x));
    for ii=1:size(x,1)
        for jj=1:size(x,2)
            z(ii,jj)=f([x(ii,jj);y(ii,jj)]);
        end
    end
    figure;
    surf(x,y,z,'edgecolor','none');
  %  shading interp
    axis off
end

end %main function

function H=compute_Hessian(x,a,theta,delta)
Diag=theta^2*a.*cos(theta*x)+2*delta;
H=diag(Diag);
end