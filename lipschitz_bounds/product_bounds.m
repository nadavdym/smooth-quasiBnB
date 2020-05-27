function [Lprod] = product_bounds(L1,L2)
%compute bounds for 0,1,2,3 order of f=f1*f2 from the individual bounds of
%each function

fix=1; %to make matlab notation consistent with "natural notation" of derivative orders
Prod=@(t,s)  L1(s+fix)*L2(t+fix);

Lprod=nan(4,1);
Lprod(0+fix) = Prod(0,0);
Lprod(1+fix)=Prod(1,0)+Prod(0,1);
Lprod(2+fix)=Prod(2,0)+Prod(0,2)+2*Prod(1,1);
Lprod(3+fix)=Prod(3,0)+Prod(0,3)+3*Prod(2,1)+3*Prod(1,2);


end

