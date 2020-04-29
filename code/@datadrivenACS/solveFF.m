function [K,obj] = solveFF(lpvdd,W,verbose)

tic
% format user specified data
format(lpvdd,'FF');
formatKFF(lpvdd);

lpvdd.internals.FF.disp = '';
lpvdd.FF.controller.Fy = tf(1);
%% 
if nargin == 2
    % verbose not specified
    verbose = 1;
elseif nargin == 1
    W = lpvdd.M.freqs;
    verbose = 1;
end


%%
% remove frequencies where very large values occur (eg at 1/s, s=1j*0)
W = reshape(W,1,length(W));
flds = fields(lpvdd.FF.constraints);
for ii = 1 : length(flds)
    if ~isempty(lpvdd.FF.constraints.(flds{ii}))
        check_Wii = (lpvdd.FF.constraints.(flds{ii})(W,1));
        W(abs(check_Wii)>1e10) = [];  
    end
end

flds = fields(lpvdd.FF.objective);
for ii = 1 : length(flds)
    if ~isempty(lpvdd.FF.objective.(flds{ii}))
        check_Wii = (lpvdd.FF.objective.(flds{ii})(W,1));
        W(abs(check_Wii)>1e10) = [];  
    end
end
%%

obj = [];
lpvdd.FF.internals.nIter = 0;
prev_obj = Inf;

lpvdd.FF.sol_obj = [];
cons = fields(lpvdd.FF.constraints);

lpvdd.FF.internals.l1 = 0;
lpvdd.FF.internals.l1_hasChanged = 0;
for i = 1 : length(cons)
   if ~isempty(lpvdd.FF.constraints.(cons{i}))
       lpvdd.FF.internals.l1  = 1;
   end
end

objs =  fields(lpvdd.FF.objective);
lpvdd.FF.internals.l2 = 1;
for i = 1 : length(objs)
   if ~isempty(lpvdd.FF.objective.(objs{i}))
       lpvdd.FF.internals.l2  = 0;
   end
end

sep ='|------|-----------|------------|----------|----------|----------|----------|------------|\n';

if verbose
    displogger(lpvdd,'FF','\n<strong>FEEDFORWARD (start: %s)</strong>\n',char(datetime(now,'ConvertFrom','datenum')));
    
    fprintf('- Solving for %d models at %d frequencies\n',lpvdd.M.nmod, length(W));
    lpvdd.FF.internals.disp = sprintf('- %d models, %d frequencies each\n',lpvdd.M.nmod, length(W));
    displogger(lpvdd,'FF','-- Controller order %d\n', lpvdd.FF.controller.order);
    if (lpvdd.Feedforward.parameters.lpv.enable) && ~(isempty(lpvdd.FF.parameters.lpv.theta))
             displogger(lpvdd,'FF','--- LPV controller -  %d scheduling parameters\n', length(theta(lpvdd,1,'FF'))-1);
     end
    sep2 ='|----------------------   COMPUTING K SATISFYING CONSTRAINTS  ---------------------------|\n';
    displogger(lpvdd,'FF',sep);
    displogger(lpvdd,'FF','| iter | objective | rel change |  slack   |  primal  |   dual   |   gap    | solve time |\n');
    displogger(lpvdd,'FF',sep2);
    
    if lpvdd.FF.internals.l2 && ~lpvdd.FF.internals.l1
        displogger(lpvdd,'FF','|-----------------------------  NO OBJECTIVE / CONSTRAINTS  -----------------------------|\n')
        K = [];
        obj = [];
        return
    end
    
end
for iter = 1 : lpvdd.Feedforward.parameters.maxIter  
    %%
     % ----- SOLVING ------ 
      t2 = tic;
      
      freqs = [];
      for j = 1 : lpvdd.M.nmod
          den = flip(lpvdd.FF.controller.K0mat.den*theta(lpvdd,j))';
          num = flip(lpvdd.FF.controller.K0mat.num*theta(lpvdd,j))';

  
          Ts = lpvdd.FF.controller.Ts;
          
          p = roots(den);
          z = roots(num);
          
          p(imag(p)<0) = [];
          z(imag(z)<0) = [];
          
          fr1 = sort(unique(imag(log(p))))/Ts;
          fr2 = sort(unique(imag(log(z))))/Ts;
          
          fr1(abs(fr1)<W(1)) = [];
          fr2(abs(fr2)<W(1)) = [];

          fr1(abs(fr1)>W(end)) = [];
          fr2(abs(fr2)>W(end)) = [];
          freqs = [freqs;fr1(:);fr2(:)];
      end
      W_ = sort(unique([W(:);freqs(:)])).';
      

      
     [Kred,obj] = lpv_solve_ff_msk9(lpvdd,W_);
     dt = toc(t2);
    
     
     
   if lpvdd.FF.internals.l1_hasChanged
        lpvdd.FF.internals.l1_hasChanged = 0;
        prev_obj = inf;
    end
    
   
    % ------------------------------
    % Check if switching from constraints <-> objective
    lpvdd.FF.sol_obj = obj;
    
    %% 
    % ----- DISPLAY AND PREMATURE EXIT ------ %
    if verbose
        displogger(lpvdd,'FF','| %4d | %4.3e | %10.2e | %4.2e | %8s | %8s | %4.2e | %4.2e s |\n',...
            iter, obj.obj, ( -obj.obj +prev_obj)/(prev_obj), max(obj.slack),obj.primal(1:7),obj.dual(1:7),...
            abs(obj.primalVal-obj.dualVal),dt);
    end
    
    
    % CAN NOT SOLVE
    if ( ~strcmpi(obj.primal,'optimal') || ~strcmpi(obj.dual,'optimal'))
        % OPTIMAL + UNKNOWN accepted
        if  ( strcmpi(obj.primal,'unknown') && strcmpi(obj.dual,'unknown'))
        else
            displogger(lpvdd,'FF',sep);
            displogger(lpvdd,'FF','--- <strong>terminating: Unable to solve at iteration %d</strong>\n',iter);
            %K = lpvdd.K0*lpvdd.Fx/lpvdd.FF.controller.Fy;
            break
        end
    end
    eigMax = 0;
    for j = 1 : lpvdd.M.nmod 
        theta_ = theta(lpvdd,j,'FF');
        den = flip(Kred.den*theta_)';
        eigMax = max(eigMax,max(abs(roots(den))));
    end
    if lpvdd.FF.parameters.lpv.enable
        for j = 1 : size(lpvdd.Feedforward.parameters.lpv.additionalPts,1)
        theta_ = lpvdd.Feedforward.parameters.lpv.theta(lpvdd.Feedforward.parameters.lpv.additionalPts(j,:));
        theta_ = theta_(:);
        den = flip(Kred.den*theta_)';
        eigMax = max(eigMax,max(abs(roots(den))));
        end
    end

    % DESTABILIZING CONTROLLER
    if eigMax >1
        displogger(lpvdd,'FF',sep);
        displogger(lpvdd,'FF','|----------------- <strong>terminating : Destabilizing controller at iteration %02d --------------|</strong>\n',iter);
        break
    end
    %%
    % ----- "VALID" CONTROLLER ------ %
    lpvdd.FF.controller.K0mat = Kred;

    
    lpvdd.FF.internals.nIter = iter;
    
    if abs( obj.obj) < lpvdd.Feedforward.parameters.exit
        if verbose
            displogger(lpvdd,'FF','|--- <strong>terminating: obj[%d]< %f</strong> (exit)\n',iter,lpvdd.Feedforward.parameters.exit);
        end
        break
    end
    
    if ( -obj.obj +prev_obj)/(obj.obj) < lpvdd.Feedforward.parameters.tol
        if verbose
            displogger(lpvdd,'FF','|--------------  <strong>terminating : |obj[%02d]-obj[%02d]|/obj[%02d]< %4.2e</strong> (tol) ----------------|\n',iter,iter-1,iter,lpvdd.Feedforward.parameters.tol);
        end
        break
    else
        if iter == lpvdd.Feedforward.parameters.maxIter
            if verbose
                displogger(lpvdd,'FF','|----------------------------- <strong>terminating : maxIter reached</strong> ----------------------------|\n');
            end
        end
    end
    prev_obj = obj.obj;
    
    lpvdd.Logs(end+1) = copy(lpvdd);
    lpvdd.Logs(end).Logs = []; % strip the logs' logs.
end % END FOR

if verbose
    t = toc;
    displogger(lpvdd,'FF','--- Finished (%s). Total run time: %f seconds.\n',char(datetime(now,'ConvertFrom','datenum')), t);
    
end

%% 
% Controller MAT OUT

lpvdd.Logs(end+1) = copy(lpvdd);
lpvdd.Logs(end).Logs = []; % strip the logs' logs.
    
K = lpvdd.FF.controller.K0mat;
lpvdd.Feedforward.sol.obj = obj;
lpvdd.Feedforward.sol.K0 = K;
end % END SOLVE
