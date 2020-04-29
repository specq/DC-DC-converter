function lpvdd = formatKFF(lpvdd)
if isa(lpvdd.Feedforward.controller.K0,'double')
    lpvdd.Feedforward.controller.K0 = tf(lpvdd.Feedforward.controller.K0,1,lpvdd.Feedforward.controller.Ts);
end



num = lpvdd.Feedforward.controller.K0.num{:};
den = lpvdd.Feedforward.controller.K0.den{:};

while ~den(1)
    % remove extra zeros (when length(num)>length(den) from removing the
    % fixed part in the numerator)
    den(1) = [];
    num = [0,num];
end

theta_ = theta(lpvdd,1,'FF');
schedParams = length(theta_);


NUM=  zeros(1+lpvdd.Feedforward.controller.order,schedParams);
DEN=  zeros(1+lpvdd.Feedforward.controller.order,schedParams);


DEN((end-length(den)+1):end,1) = flip(den');
NUM((end-length(num)+1):end,1) = flip(num');


lpvdd.FF.controller.K0mat = [];
lpvdd.FF.controller.K0mat.num = NUM;
lpvdd.FF.controller.K0mat.den = DEN;


lpvdd.FF.controller.K0mat.Ts = lpvdd.Feedforward.controller.Ts;
end