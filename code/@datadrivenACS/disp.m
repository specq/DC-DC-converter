function disp(lpvdd)
fprintf('LPV Data-driven optimization problem \n')
% fprintf('Accessible properties \n')
% props = properties(lpvdd);
% 
% props_ = pad(props) ;
% blk = repmat(' ',1,length(props_{1}));

% for ii = 1 : length(props) 
%     f = fields(lpvdd.(props{ii}));
%     for jj = 1 : length(f)
%         if jj == 1
%            fprintf('<strong>%s</strong>',props_{ii});
%            fprintf(' ')
%            fprintf(f{jj}) 
%            
%         else
%            fprintf(blk) 
%            fprintf(' ')
%            fprintf(f{jj}) 
%         end
%         fprintf('\n')
%     end
%     fprintf('\n')
% end


if numel(lpvdd)==1
    
    fprintf(lpvdd.internals.FB.disp);
    fprintf(lpvdd.internals.FF.disp);
    fprintf('\n')
else
    fprintf('\n')
    for jj = 1 : numel(lpvdd)
        str_disp = strcat('%',sprintf('%dd [datadriven]',ceil(log10(numel(lpvdd)))),'\n'); 
        fprintf(str_disp,jj);
    end
    
end

end % END DISPLAY


