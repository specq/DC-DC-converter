function [z,i] = exp_sol_int(x) %#codegen
%Evaluate function "primal" with tiebreak "first-region"
% 
%  [value, region] = exp_sol_int(x)
%
%See "help PolyUnion/toMatlab" for more information.
x=x(:);xh=[x;-1];
if numel(x)~=2,error('The input vector must have 2 elements.');end
nx=2;nz=1;
H=[-1 0 5000;0.705168472560137 -0.70903979106056 109.50190555189;-0.0216520381544553 -0.999765567142497 4802.33619378047;0.616356115007109 0.787467548216016 -6196.6282392853;-0.705168472560137 0.70903979106056 -109.501905551889;-0.00561180785388593 -0.999984253682332 4729.98213546524;1 0 10000;0.00372659958448609 0.99999305620366 -4487.65309227231;-1 0 5000;0.705168472560137 -0.70903979106056 109.50190555189;-0.0235895996715379 -0.999721726675647 5081.62361750681;0.0216520381544553 0.999765567142497 -4802.33619378047;-0.616356115007109 -0.787467548216016 6196.6282392853;0.616356115007109 0.787467548216016 3163.83737897044;-1 0 5000;-0.00372659958448609 -0.99999305620366 4487.65309227231;1 0 10000;0.00372659958448609 0.99999305620366 2243.82654613615;-0.705168472560137 0.70903979106056 -109.50190555189;0 -1 5000;-0.00561180785388596 -0.999984253682332 4995.9010403669;1 0 10000;0.00561180785388593 0.999984253682332 -4729.98213546524;-1 0 5000;0 -1 5000;0.705168472560137 -0.70903979106056 109.501905551891;0.0235895996715379 0.999721726675647 -5081.62361750681;0.966654626584156 0.256083644349744 9963.42247151394;1 0 10000;0.0216520381544553 0.999765567142497 2402.87835201801;-0.705168472560137 0.70903979106056 -20.8545959395862;-0.616356115007109 -0.787467548216016 -3163.83737897044;-0.00372659958448609 -0.99999305620366 -2243.82654613615;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388588 0.999984253682332 2364.99106773262;-0.705168472560137 0.70903979106056 -109.50190555189;0 -1 5000;0.00561180785388596 0.999984253682332 -4995.9010403669;1 0 10000;0.0403705975114638 0.999184775132491 2708.67351350464;-0.4389862378058 0.898493785742066 -2299.45327720639;-0.966654626584156 -0.256083644349743 -9963.42247151394;-0.0216520381544553 -0.999765567142497 -2402.87835201801;0.438986237805801 -0.898493785742066 2299.45327720639;0.0235895996715379 0.999721726675647 2542.72866264449;-0.705168472560137 0.70903979106056 -20.8545959395872;-0.00561180785388588 -0.999984253682332 -2364.99106773262;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388591 0.999984253682332 2497.95052018345;-0.0403705975114639 -0.999184775132491 -2708.67351350464;1 0 10000;-0.438986237805801 0.898493785742066 -2299.45327720639;-0.0235895996715379 -0.999721726675647 -2542.72866264449;0.438986237805801 -0.898493785742066 2299.45327720639;1 0 10000;0.0235895996715393 0.999721726675647 2680.47697485512;-0.705168472560137 0.70903979106056 -20.8545959395875;-0.00561180785388591 -0.999984253682332 -2497.95052018345;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388564 0.999984253682332 2633.30119195576;-0.0235895996715393 -0.999721726675647 -2680.47697485512;1 0 10000;0.0235895996715387 0.999721726675647 2820.7026319975;-0.705168472560137 0.70903979106056 -20.8545959395871;-0.00561180785388564 -0.999984253682332 -2633.30119195576;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.005611807853886 0.999984253682332 2771.08608811415;-0.0235895996715387 -0.999721726675647 -2820.70263199751;1 0 10000;0.0235895996715385 0.999721726675647 2963.4501880688;-0.705168472560137 0.70903979106056 -20.8545959395872;-0.005611807853886 -0.999984253682332 -2771.08608811415;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388597 0.999984253682332 2911.34898715107;-0.0235895996715385 -0.999721726675647 -2963.4501880688;1 0 10000;0.0235895996715389 0.999721726675647 3108.76499835093;-0.705168472560137 0.70903979106056 -20.8545959395863;-0.00561180785388597 -0.999984253682332 -2911.34898715107;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388606 0.999984253682332 3054.13445489661;-0.0235895996715389 -0.999721726675647 -3108.76499835093;1 0 10000;0.0235895996715388 0.999721726675647 3256.69323382124;-0.705168472560137 0.70903979106056 -20.8545959395867;-0.00561180785388606 -0.999984253682332 -3054.13445489661;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388599 0.999984253682332 3199.48785867831;-0.0235895996715388 -0.999721726675647 -3256.69323382124;1 0 10000;0.0235895996715387 0.999721726675647 3407.28189582256;-0.705168472560137 0.70903979106056 -20.8545959395875;-0.00561180785388599 -0.999984253682332 -3199.48785867831;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388597 0.999984253682332 3347.45538173592;-0.0235895996715387 -0.999721726675647 -3407.28189582256;1 0 10000;0.0235895996715388 0.999721726675647 3560.57883099691;-0.705168472560137 0.70903979106056 -20.8545959395869;-0.00561180785388597 -0.999984253682332 -3347.45538173592;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.005611807853886 0.999984253682332 3498.08403789512;-0.0235895996715388 -0.999721726675647 -3560.57883099691;0.0235895996715388 0.999721726675647 3716.63274648785;1 0 10000;-0.705168472560137 0.70903979106056 -20.8545959395869;-0.00561180785388599 -0.999984253682332 -3498.08403789512;0.705168472560137 -0.70903979106056 20.8545959395869;-1 0 5000;0.00561180785388607 0.999984253682332 3651.4216865053];
ni=[1;5;9;13;19;24;28;33;37;40;44;48;52;55;60;64;68;72;76;80;84;88;92;96;100;104;108;112;116;120];
fF=[0 0;-0.0652885851869618 0.0656470142965094;0 0;-0.0658467367056002 -0.0841269633724433;-0.0652885851869618 0.0656470142965094;0 0;0 0;-0.0652885851869618 0.0656470142965094;-0.0652885851869618 0.0656470142965094;1.38777878078145e-17 -1.38777878078145e-17;0 -1.38777878078145e-17;-0.0652885851869618 0.0656470142965094;0 0;1.38777878078145e-17 0;-0.0652885851869618 0.0656470142965094;0 -1.38777878078145e-17;-0.0652885851869618 0.0656470142965094;0 0;-0.0652885851869618 0.0656470142965094;0 -1.38777878078145e-17;-0.0652885851869618 0.0656470142965094;0 0;-0.0652885851869618 0.0656470142965094;1.38777878078145e-17 1.38777878078145e-17;-0.0652885851869618 0.0656470142965094;0 0;-0.0652885851869618 0.0656470142965094;0 0;-0.0652885851869618 0.0656470142965094];
fg=[662;672.138321219614;662;0;672.138321219614;662;-338;-336.069160609807;672.138321219614;-338;-338;-336.069160609807;-338;-338;-336.069160609807;-338;-336.069160609807;-338;-336.069160609807;-338;-336.069160609807;-338;-336.069160609807;-338;-336.069160609807;-338;-336.069160609807;-338;-336.069160609807];
for i=1:29,
if all(H(ni(i):ni(i+1)-1,:)*xh<=1e-08);
z=fF((i-1)*nz+1:i*nz,:)*x+fg((i-1)*nz+1:i*nz);
return
end
end
i=0;z=NaN(1,1);
end