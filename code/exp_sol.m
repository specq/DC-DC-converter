function [z,i] = exp_sol(x) %#codegen
%Evaluate function "primal" with tiebreak "obj"
% 
%  [value, region] = exp_sol(x)
%
%See "help PolyUnion/toMatlab" for more information.
x=x(:);xh=[x;-1];
if numel(x)~=2,error('The input vector must have 2 elements.');end
nx=2;nz=10;
H=[0.999456236344815 -0.0329731956512758 0.543916341112487;-1 0 0.05;0.995484759334458 0.094921514593991 0.119230377256846;0 -1 5;1 0 0.45;-0.961311217383224 0.275464595425953 -1.3292574162606;0 -1 5;-0.999456236344815 0.0329731956512758 -0.543916341112487;0.961311217383224 -0.275464595425953 1.3292574162606;1 0 0.45;-1 0 0.05;0.995535067130876 0.0943924261406664 0.267817809440488;-0.995484759334458 -0.094921514593991 -0.119230377256846;1 0 0.45;0 1 5;-1 0 0.05;0.995535067130877 0.0943924261406665 0.452285636184273;-0.995535067130877 -0.0943924261406663 -0.267817809440489;1 0 0.45;0 1 5;0.995535067130877 0.0943924261406665 0.642718821678027;-0.995535067130877 -0.0943924261406665 -0.452285636184273];
ni=[1;5;9;14;19;23];
fF=[-0.618391731398726 0.0204014451129071;0.344089771253342 0.032716128387993;-0.00583888658907661 -0.000706748082138386;-0.000539422804310519 -4.81947433588897e-05;2.21846495137012e-05 2.33631711617222e-06;5.7969151501247e-07 4.4563131795794e-08;-5.4292286222335e-08 -5.38119579104879e-09;-1.48470125083122e-12 2.45303465040703e-11;1.03945768659131e-10 9.79898454045447e-12;-2.11849982001411e-12 -2.46937886805298e-13;0 0;0.0885375109761127 0.0411470875143911;0.4242321246834 -0.0148952788881406;-0.0088275655513782 0.000225240493169548;-0.000631819697509362 2.39126628399644e-05;2.97906999165432e-05 -9.1914119078329e-07;6.0130819901083e-07 -2.70101999273742e-08;-6.92914424171853e-08 2.3104846993971e-09;2.63257027288688e-10 4.5431297612808e-12;1.27251098547276e-10 -4.51498560760655e-12;-0.20781815132297 0.0595504784935119;0.174418992892269 0.0165376714971594;0.279701468851234 0.0265201107541752;-0.00604223318656705 -0.000572899005387505;-0.0004120334884537 -3.90672733793715e-05;1.99739810697874e-05 1.89384843903018e-06;3.80985588077642e-07 3.61234427335988e-08;-4.60056995277291e-08 -4.36206593415478e-09;2.09718464816433e-10 1.98846286658672e-11;8.37748481696821e-11 7.94317389640753e-12;-0.20781815132297 0.0595504784935119;0.509117585916495 0.0482723770484917;0.141386259225963 0.0134056473467623;0.226729576436146 0.0214975398700245;-0.00489791125068206 -0.000464399237392998;-0.000333999598565671 -3.16684298519582e-05;1.61911636942236e-05 1.53517768846506e-06;3.08831774875618e-07 2.92821232053542e-08;-3.72928065095479e-08 -3.53594624413756e-09;1.70000458155073e-10 1.61187209413161e-11;-0.20781815132297 0.0595504784935119;0.509117585916495 0.0482723770484917;0.412697205649767 0.0391302042378321;0.114609504196933 0.0108667886417156;0.18378988513019 0.0174261798808865;-0.00397030930101955 -0.000376447942242897;-0.000270744332604655 -2.56708329642585e-05;1.3124763704353e-05 1.24443460554347e-06;2.50342973684425e-07 2.3736462372026e-08;-3.0230024372635e-08 -2.86628309628667e-09];
fg=[4.9960036108132e-16;1.11022302462516e-16;2.22044604925031e-16;5.55111512312578e-16;4.44089209850063e-16;0;4.44089209850063e-16;1.11022302462516e-16;4.44089209850063e-16;4.44089209850063e-16;-0.336536364160089;0.139074674126163;-0.234049918709165;0.00451050892842508;0.000355917186334986;-1.58969890151628e-05;-3.56785823774075e-07;3.77084449221243e-08;-8.66988703052129e-11;-7.04041269727895e-11;-0.0491748792586778;0.0203216681358118;-0.034199503289195;0.000659078053880546;5.20068156910636e-05;-2.32287680801413e-06;-5.21337408798672e-08;5.50997880655046e-09;-1.26683108447878e-11;-1.02873265461767e-11;-0.0491748792586778;-0.0697185996218056;0.00300991070051071;-0.0619609523307026;0.00125879299064757;9.29027278218086e-05;-4.30537142936949e-06;-8.99480302463118e-08;1.00762218568917e-08;-3.3483438244275e-11;-0.0491748792586778;-0.0697185996218055;-0.12025048268844;-0.0110232205761739;-0.0844647368416895;0.00174492956388006;0.00012605347555894;-5.91240684399974e-06;-1.20600774033974e-07;1.37776765551934e-08];
tH=[135.682534103826 3.33425643971433;3.33425643971433 1.31775701387808;273.33090557106 -1.20691957022027;-1.20691957022027 1.46757556482341;196.359761655315 9.11995459987719;9.11995459987719 1.86943520603192;236.682637919419 12.9431992501214;12.9431992501214 2.23193909984848;263.178518785315 15.4554266544999;15.4554266544999 2.47013788113278];
tF=[-3.57523017875034e-15 -3.43343955629704e-16;-149.819863733824 4.9427273548359;-14.5347654276499 -1.38591970969329;-36.2300019852435 -3.44297031692994;-60.3049074895399 -5.72565107292491];
tg=[3.73286533643419e-28;40.7669036145563;0.870422951747174;3.78863795088998;9.25742268364088];
tb=[];
for i=1:5,
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
i=0;z=NaN(10,1);
end
