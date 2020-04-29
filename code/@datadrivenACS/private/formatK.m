function lpvdd = formatK(lpvdd)

% if isa(lpvdd.Feedback.controller.K0,'tf')
%     % if K0 is a TF, the sampling time and poles close to the unit circle can be retrived
%     if lpvdd.Feedback.controller.K0.Ts
%         K = tf(lpvdd.Feedback.controller.K0);
%         setKinit(lpvdd,K);
%     end
% end

if ~lpvdd.Feedback.parameters.lpv.K0mat
    % lpvdd.Feedback.parameters.K0mat == 1 -> setCustomK0mat has been used
    % to set the inital controller, obtain by solving the optimization
    % problem earlier.
    
    if isa(lpvdd.Feedback.controller.K0,'double')
        % convert to transfer function
        
        if isempty(lpvdd.Feedback.controller.Ts)
            error('Missing sampling time')
        end
        lpvdd.Feedback.controller.K0 = tf(squeeze(lpvdd.Feedback.controller.K0)',1,lpvdd.Feedback.controller.Ts);
    end
    
    if isa(lpvdd.Feedback.controller.Fy,'double')
        % convert to transfer function
        lpvdd.Feedback.controller.Fy = tf(squeeze(lpvdd.Feedback.controller.Fy)',1,lpvdd.Feedback.controller.Ts);
        lpvdd.FB.controller.Fy = lpvdd.Feedback.controller.Fy;
    end
    
    num = lpvdd.Feedback.controller.K0.num{:};
    den = lpvdd.Feedback.controller.K0.den{:};
    
    while ~den(1)
        % remove extra zeros (when length(num)>length(den) from removing the
        % fixed part in the numerator)
        den(1) = [];
        % num = [0,num];
    end
    
    % get size of theta
    theta_ = theta(lpvdd);
    schedParams = length(theta_);
    
    sx = 1+lpvdd.Feedback.controller.order;
    sx = max(sx,size(num,2));
    
    sy = 1+lpvdd.Feedback.controller.order-order(1/lpvdd.Feedback.controller.Fy);
    sy = max(sy,size(den,2));
    
    NUM=  zeros(sx,schedParams);
    DEN=  zeros(sy,schedParams);
    
    
    DEN((end-length(den)+1):end,1) = flip(den');
    NUM((end-length(num)+1):end,1) = flip(num');
    
    lpvdd.FB.controller.K0mat = [];
    lpvdd.FB.controller.K0mat.num = NUM;
    lpvdd.FB.controller.K0mat.den = DEN;
    lpvdd.FB.controller.K0mat.Ts = lpvdd.Feedback.controller.Ts;
end
end
