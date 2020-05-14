function [z,i] = exp_sol_2(x) %#codegen
%Evaluate function "primal" with tiebreak "first-region"
% 
%  [value, region] = exp_sol_2(x)
%
%See "help PolyUnion/toMatlab" for more information.
x=x(:);xh=[x;-1];
if numel(x)~=2,error('The input vector must have 2 elements.');end
nx=2;nz=1;
H=[1 0 0.222222222222222;0 1 5;-1 0 0.277777777777778;0 -1 5];
ni=[1;5];
fF=[-0.529982214636897 0.0250510149313122];
fg=0;
for i=1:1,
if all(H(ni(i):ni(i+1)-1,:)*xh<=1e-08);
z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);
return
end
end
i=0;z=NaN(1,1);
end
