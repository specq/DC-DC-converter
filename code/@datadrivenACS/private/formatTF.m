function dd = formatTF(dd,in1,in2,in3)

flds = fields(dd.(in1).(in3));
for ii = 1 : length(flds)
    arg2 = dd.(in1).(in3).(flds{ii});
    if ~isempty(arg2)
                arg2 = arg2(:,:);
                len_ = length(arg2);
                switch class(arg2)
                    case 'function_handle'
                        f_handle_arg2 = arg2;
                    case {'tf';'frd';'zpk';'ss';'idpoly'}
                        if len_ > 1
                            f_handle_arg2 = @(w,m) freqresp(arg2(:,:,m),w);
                        else
                            f_handle_arg2 = @(w,~) freqresp(arg2,w); 
                        end
                    case 'double'
                        f_handle_arg2 = @(w,~) arg2;
                end
                dd.(in2).(in3).(flds{ii}) = f_handle_arg2;
    end              
end

end
