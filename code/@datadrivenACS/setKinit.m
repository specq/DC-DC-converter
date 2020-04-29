function dd = setKinit(dd,Kinit)

Kinit = tf(Kinit); % in case if Kinit is other than TF

Ts_ = Kinit.Ts;
Fy_ = tf(1);
if Ts_ > 0
    g = eig(Kinit);
    if ~isempty(g) % found some poles
        z = tf('z',Ts_);
        fy = g((abs(abs(g)-1)<1e-4));
        for  ii = 1:length(fy)
            Fy_ = Fy_*(z-fy(ii));
        end
    end
    
        dd.Feedback.controller.Fy = Fy_;
   
    dd.Feedback.controller.K0 = minreal(Kinit*Fy_,1e-6);
    dd.Feedback.controller.Ts = Kinit.Ts;
    dd.Feedforward.controller.Ts = Kinit.Ts;
end % END ADD POLES FY DISCRETE
end

