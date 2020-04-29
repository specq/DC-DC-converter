function [Klpv,KlpvFF] = computeK(lpvdd,arg2,theta_in)
%{
computeK(lpvdd) to sample the LPV controller at ALL the different
SamplingGrid points

computeK(lpvdd,theta,'FB'/'FF') to get the LPV controller sampled at theta.

[KFB,KFF] = computeK(P) to get FB and FF controller

%}
if nargin < 3
    theta_in = [];
end

if nargin < 2
    arg2 = 'FB';
end

if strcmpi(arg2,'Feedback')
    arg2 = 'FB';
end

if strcmpi(arg2,'Feedforward')
    arg2 = 'FF';
end

if ~isempty(theta_in)
    % nummber of controller to get. 
    m = 1;
else
    m = lpvdd.M.nmod;
end

for j = 1 : m
    theta_ = theta_in;
    if isempty(lpvdd.(arg2).parameters.lpv.theta) ||  lpvdd.(arg2).parameters.lpv.enable==0
        % Not a LPV controller
        theta_ = 1;
    else
        if ~isempty(theta_)
            theta_ = lpvdd.(arg2).parameters.lpv.theta(theta_);
        else
            theta_ = theta(lpvdd,j,arg2);
        end
    end
    theta_ = theta_(:);
    
    num = flip(lpvdd.(arg2).controller.K0mat.num*theta_)';
    den = flip(lpvdd.(arg2).controller.K0mat.den*theta_)';
    
    if  j == 1
        Klpv= tf(num,den,lpvdd.(arg2).controller.Ts)/lpvdd.(arg2).controller.Fy;
    else
        Klpv(:,:,j) = tf(num,den,lpvdd.(arg2).controller.Ts)/lpvdd.(arg2).controller.Fy;
    end
end


if size(lpvdd.(arg2).controller.K0mat.num,2) == 1 || nargin == 2
    Klpv = Klpv(:,:,1);
end

if nargout == 2
    % in case we want to fetch both FB and FF
    KlpvFF = computeK(lpvdd,'FF',theta_in);
end

end
