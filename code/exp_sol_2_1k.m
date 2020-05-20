function [z,i] = exp_sol_2_1k(x) %#codegen
%Evaluate function "primal" with tiebreak "first-region"
% 
%  [value, region] = exp_sol_2_1k(x)
%
%See "help PolyUnion/toMatlab" for more information.
x=x(:);xh=[x;-1];
if numel(x)~=2,error('The input vector must have 2 elements.');end
nx=2;nz=1;
H=[-1 0 0.277777777777778;0 -1 5;0.99761733903292 0.0689901794524095 -0.555411963936045;1 0 0.122222222222222;-0.99761733903292 -0.0689901794524095 0.555411963936045;0 1 5;-1 0 0.277777777777778;0 -1 5];
ni=[1;4;9];
fF=[-0.314257478347589 0.0388866896794344;-0.520895186167106 0.0245966688605073];
fg=[0.115043163979625;3.33066907387547e-16];
for i=1:2,
if all(H(ni(i):ni(i+1)-1,:)*xh<=1e-08);
z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);
return
end
end
i=0;z=NaN(1,1);
end
