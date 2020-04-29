function lpvdd = displogger(lpvdd,arg2,varargin)
% Some operations have to be done after the different
% properties have been specified. Example scheduling:
% we have to both know the initial controller, controller order
% and scheduling order to construct the correct intinial matrix
% coeffeicient gains

newline = sprintf(varargin{:});

fprintf(newline);
lpvdd.internals.(arg2).disp = strcat(lpvdd.internals.(arg2).disp,'\n',newline);

end
