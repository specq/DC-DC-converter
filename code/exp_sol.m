function [z,i] = exp_sol(x) %#codegen
%Evaluate function "primal" with tiebreak "obj"
% 
%  [value, region] = exp_sol(x)
%
%See "help PolyUnion/toMatlab" for more information.
x=x(:);xh=[x;-1];
if numel(x)~=2,error('The input vector must have 2 elements.');end
nx=2;nz=1;
H=[1 0 0.222222222222222;0 1 5;-1 0 0.277777777777778;0 -1 5];
ni=[1;5];
fF=[-0.529982214636897 0.0250510149313122];
fg=0;
tH=[114.417357205858 0.936120031293399;0.936120031293399 1.06440493495907];
tF=[1.4716723152044e-16 1.30456018054158e-17];
tg=1.74783280394215e-28;
tb=[];
for i=1:1,
if all(H(ni(i):ni(i+1)-1,:)*xh<=1e-08);
tv=tF(i,:)*x+tg(i);
tv=tv+x'*tH((i-1)*nx+1:i*nx,:)*x;
tb=[tb;i,tv];
end
end
if ~isempty(tb)
[~,j]=min(tb(:,end));
i=tb(j,1);
z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);
return
end
i=0;z=NaN(1,1);
end
