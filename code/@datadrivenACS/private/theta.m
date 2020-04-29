function  theta_out = theta(dd,mod,arg3)
if nargin == 1
    mod = 1;
end
if nargin < 3
    arg3 = 'FB';
end

if strcmpi(arg3,'fb')
    arg3 = 'Feedback';
end

if strcmpi(arg3,'ff')
    arg3 = 'Feedforward';
end

%[~,sz] = size(dd.M.SamplingGrid);

if dd.(arg3).parameters.lpv.enable && ~isempty(dd.(arg3).parameters.lpv.theta)
        theta_out = dd.(arg3).parameters.lpv.theta(dd.M.SamplingGrid(mod,:));
else
    theta_out = 1;
end
theta_out = theta_out(:);
end % END ADD