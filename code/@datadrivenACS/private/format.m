function dd = format(dd,arg2)
dd = formatModel(dd);

if nargin == 1
    arg2 = 'FB';
end

if strcmpi(arg2,'FB')
    if dd.Feedback.parameters.lpv.K0mat
        % special case, whem K0mat (LPV) has been specified
       tmp = dd.FB.controller.K0mat;  
    end
    dd.FB = dd.Feedback;
    
    if dd.Feedback.parameters.lpv.K0mat
       dd.FB.controller.K0mat = tmp;  
    end
    dd = formatTF(dd,'Feedback','FB','constraints');
    dd = formatTF(dd,'Feedback','FB','objective');
    
elseif strcmpi(arg2,'FF')
    dd.FF = dd.Feedforward;
    dd = formatTF(dd,'Feedforward','FF','constraints');
    dd = formatTF(dd,'Feedforward','FF','objective');
end

end

