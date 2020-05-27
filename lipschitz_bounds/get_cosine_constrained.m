function [c0,f,g,H,L1,L2,L3,gt]=get_cosine_constrained()
d=2;
alpha=10*ones(d,1);
k=1;
theta=2*k*pi;
seed=1;
rng(seed);
%W=diag([1 -2]);
W=randn(2);
W=W+W';

gt=inf(d,1);

f=@(x) alpha'*(1-cos(theta*x))+x'*W*x;

c0.x=zeros(d,1);
c0.h=ones(d,1);

L1=theta*norm(alpha)+2*norm(W)*norm(c0.h);
L2=theta^2*norm(alpha)+2*norm(W);
L3=norm(alpha)*theta^3;

g=@(x) theta*alpha.*sin(theta*x)+2*W*x;
H=@(x) compute_Hessian(x,alpha,theta,W);


%plot function
doDisplay=true;
if doDisplay && d==2
    [x,y]=meshgrid(-1:0.05:1);
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

function H=compute_Hessian(x,a,theta,W)
Diag=theta^2*a.*cos(theta*x);
H=diag(Diag)+2*W;
end