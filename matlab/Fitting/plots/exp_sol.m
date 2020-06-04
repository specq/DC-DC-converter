function [z,i] = exp_sol(x) %#codegen
%Evaluate function "primal" with tiebreak "first-region"
% 
%  [value, region] = exp_sol(x)
%
%See "help PolyUnion/toMatlab" for more information.
x=x(:);xh=[x;-1];
if numel(x)~=2,error('The input vector must have 2 elements.');end
nx=2;nz=1;
H=[-1 0 0.0499999999999999;-0.948234075012487 0.317572257896078 0.246980675535447;0.948234075012487 -0.317572257896078 0.740942026606345;0.997561940462708 0.0697866386944969 0.064337502802721;-0.997561940462708 -0.0697866386944969 0.126088277332991;1 0 0.15;-1 0 0.0499999999999999;0.489392926574889 0.872063394151174 1.53393273358675;0.999949453333745 -0.0100543909583595 0.000274482125110233;0.948234075012487 -0.317572257896078 -0.246980675535447;0 -1 5;-0.999949453333746 0.0100543909583595 -0.0530033380831395;-0.948234075012487 0.317572257896078 -0.740942026606345;1 0 0.15;0.920728946121259 0.390202777763609 0.698887964698173;0.973383914133187 0.229180618087036 0.447209536073511;0.987655721066322 0.15664027785012 0.353985894428786;-0.999949453333746 0.0100543909583595 -0.000274482125110234;-0.997561940462708 -0.0697866386944969 -0.064337502802721;1 0 0.15;0.997561940462708 0.0697866386944969 -0.126088277332991;-1 0 0.0499999999999999;0 -1 5;0.999949453333746 -0.0100543909583595 0.0530033380831395];
ni=[1;7;11;15;21;25];
fF=[-5.23783609496595 -0.366377373128297;-6.52665675275009 0.0655814182142277;-6.52842051937485 0.0656460197581129;0.000323259966063604 1.31119941162772e-05;0.000442908974944416 1.27559817533814e-05];
fg=[5.08647857246904e-05;-0.335745174156225;1.00811480440142;-0.337836112387269;0.662076279017365];
for i=1:5,
if all(H(ni(i):ni(i+1)-1,:)*xh<=1e-08);
z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);
return
end
end
i=0;z=NaN(1,1);
end
