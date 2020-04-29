function [K,obj] = solveFB(lpvdd,W,verbose)

tic
% format user specified data
format(lpvdd);
formatK(lpvdd);
% Trim frequencies in model
lpvdd.M.freqs = lpvdd.M.freqs(lpvdd.M.freqs<=pi/lpvdd.Feedback.controller.Ts);
lpvdd.M.freqs(end) = pi/lpvdd.Feedback.controller.Ts;


lpvdd.internals.FB.disp = '';
%% 
if nargin == 3
    if isempty(W)
        W = lpvdd.M.freqs;
    end
end
if nargin == 2
    % verbose not specified
    verbose = 1;
    
elseif nargin == 1
    W = lpvdd.M.freqs;
    verbose = 1;
end


%%
% remove very large values (eg at 1/s, s=0)
W = reshape(W,1,length(W));
flds = fields(lpvdd.FB.constraints);
for ii = 1 : length(flds)
    if ~isempty(lpvdd.FB.constraints.(flds{ii}))
        check_Wii = (lpvdd.FB.constraints.(flds{ii})(W,1));
        W(abs(check_Wii)>1e10) = [];  
    end
end

flds = fields(lpvdd.FB.objective);
for ii = 1 : length(flds)
    if ~isempty(lpvdd.FB.objective.(flds{ii}))
        check_Wii = (lpvdd.FB.objective.(flds{ii})(W,1));
        W(abs(check_Wii)>1e10) = [];  
    end
end
%%


dd.flag.largeObj = 0;

%%
obj = [];
lpvdd.FB.internals.nIter = 0;
% ------  check(lpvdd,W); % check if correctly set-up
%lpvdd.K0 = lpvdd.user_data.K; % get starting-point
prev_obj = Inf;

lpvdd.FB.sol_obj = [];
cons = fields(lpvdd.FB.constraints);

lpvdd.Feedback.sol.obj = [];
lpvdd.FB.internals.l1 = ~isempty(lpvdd.FB.parameters.c2);
lpvdd.FB.internals.l1_hasChanged = 0;
for i = 1 : length(cons)
   if ~isempty(lpvdd.FB.constraints.(cons{i}))
       lpvdd.FB.internals.l1  = 1;
   end
end

objs =  fields(lpvdd.FB.objective);
lpvdd.FB.internals.l2 = 1;
for i = 1 : length(objs)
   if ~isempty(lpvdd.FB.objective.(objs{i}))
       lpvdd.FB.internals.l2  = 0;
   end
end

sep ='|------|-----------|------------|----------|----------|----------|----------|------------|\n';

if verbose
    displogger(lpvdd,'FB','<strong>FEEDBACK (start: %s)</strong>\n',char(datetime(now,'ConvertFrom','datenum')));
    
    fprintf('- Solving for %d models at %d frequencies\n',lpvdd.M.nmod, length(W));
    lpvdd.FB.internals.disp = sprintf('- %d models, %d frequencies each\n',lpvdd.M.nmod, length(W));
    displogger(lpvdd,'FB','-- Controller order %d, numerator fixed part order %d\n', lpvdd.FB.controller.order,order(1/lpvdd.FB.controller.Fy));
    if (lpvdd.Feedback.parameters.lpv.enable) && ~(isempty(lpvdd.FB.parameters.lpv.theta))
             displogger(lpvdd,'FB','--- LPV controller -  %d scheduling parameters\n', length(theta(lpvdd,1))-1);
     end
    sep2 ='|----------------------   COMPUTING K SATISFYING CONSTRAINTS  ---------------------------|\n';
    displogger(lpvdd,'FB',sep);
    displogger(lpvdd,'FB','| iter | objective | rel change |  slack   |  primal  |   dual   |   gap    | solve time |\n');
    displogger(lpvdd,'FB',sep2);
    
    if lpvdd.FB.internals.l2 && ~lpvdd.FB.internals.l1
        displogger(lpvdd,'FB','|-----------------------------  NO OBJECTIVE / CONSTRAINTS  -----------------------------|\n')
        K = [];
        obj = [];
        return
    end
    
end
for iter = 1 : lpvdd.Feedback.parameters.maxIter  
    %%
     % ----- SOLVING ------ 
      t2 = tic;
      
      % Get frequencies for poles/zeros 
      freqs = [];
      for j = 1 : lpvdd.M.nmod
          den = flip(lpvdd.FB.controller.K0mat.den*theta(lpvdd,j))';
          num = flip(lpvdd.FB.controller.K0mat.num*theta(lpvdd,j))';

  
          Ts = lpvdd.Feedback.controller.Ts;
          
          p = roots(den);
          z = roots(num);
          
          p(imag(p)<0) = [];
          z(imag(z)<0) = [];
          
          fr1 = sort(unique(imag(log(p))))/Ts;
          fr2 = sort(unique(imag(log(z))))/Ts;
          
          fr1(abs(fr1)<W(1)+eps) = [];
          fr2(abs(fr2)<W(1)+eps) = [];

          fr1(abs(fr1)>W(end)) = [];
          fr2(abs(fr2)>W(end)) = [];
          freqs = [freqs;fr1(:);fr2(:)];
      end
      W_ = (([W(:);freqs(:)])).';


     [Kred,obj] = lpv_solve_msk9(lpvdd,W_,length(W(:)),verbose);
     dt = toc(t2);
    
     
     
   if lpvdd.FB.internals.l1_hasChanged
        lpvdd.FB.internals.l1_hasChanged = 0;
        prev_obj = inf;
    end
    
   
    % ------------------------------
    % Check if switching from constraints <-> objective
    lpvdd.FB.sol_obj = obj;
    
    %% 
    % ----- DISPLAY AND PREMATURE EXIT ------ %
    if verbose
        displogger(lpvdd,'FB','| %4d | %4.3e | %10.2e | %4.2e | %8s | %8s | %4.2e | %4.2e s |\n',...
            iter, abs(obj.obj), ( -obj.obj +prev_obj)/(prev_obj), max(obj.slack),obj.primal(1:7),obj.dual(1:7),...
            abs(obj.primalVal-obj.dualVal),dt);
    end
    
    
    % CAN NOT SOLVE
    if ( ~strcmpi(obj.primal,'optimal') || ~strcmpi(obj.dual,'optimal'))
        % OPTIMAL + UNKNOWN accepted
        if  ( strcmpi(obj.primal,'unknown') && strcmpi(obj.dual,'unknown'))
        else
            displogger(lpvdd,'FB',sep);
            displogger(lpvdd,'FB','--- <strong>terminating: Unable to solve at iteration %d</strong>\n',iter);
            %K = lpvdd.K0*lpvdd.Fx/lpvdd.FB.controller.Fy;
            break
        end
    end
    eigMax = 0;
    for j = 1 : lpvdd.M.nmod 
        theta_ = theta(lpvdd,j);
        den = flip(Kred.den*theta_)';
        eigMax = max(eigMax,max(abs(roots(den))));
    end
     if lpvdd.FB.parameters.lpv.enable
        for j = 1 : size(lpvdd.Feedback.parameters.lpv.additionalPts,1)
        theta_ = lpvdd.Feedback.parameters.lpv.theta(lpvdd.Feedback.parameters.lpv.additionalPts(j,:));
        theta_ = theta_(:);
        den = flip(Kred.den*theta_)';
        eigMax = max(eigMax,max(abs(roots(den))));
        end
    end

    % DESTABILIZING CONTROLLER
    if eigMax >1
        displogger(lpvdd,'FB',sep);
        displogger(lpvdd,'FB','|----------------- <strong>terminating : Destabilizing controller at iteration %02d --------------|</strong>\n',iter);
        break
    end
    %%
    % ----- "VALID" CONTROLLER ------ %
    lpvdd.FB.controller.K0mat = Kred;
    lpvdd.Logs(end+1) = copy(lpvdd);
    lpvdd.Logs(end).Logs = []; % strip the logs' logs.
    
    lpvdd.FB.internals.nIter = iter;
    if ( abs( obj.obj) < lpvdd.Feedback.parameters.exit && ~lpvdd.FB.internals.l1)
        if verbose
            displogger(lpvdd,'FB','|--- <strong>terminating: obj[%d]< %f</strong> (exit)\n',iter,lpvdd.Feedback.parameters.exit);
        end
        break
    end
    
    if (abs( -obj.obj +prev_obj)/(obj.obj) < lpvdd.Feedback.parameters.tol && (~(iter==1) && abs( obj.obj)>1e-8) ) 
        if verbose
            displogger(lpvdd,'FB','|--------------  <strong>terminating : |obj[%02d]-obj[%02d]|/obj[%02d]< %4.2e</strong> (tol) ----------------|\n',iter,iter-1,iter,lpvdd.Feedback.parameters.tol);
        end
        break
    else
        if iter == lpvdd.Feedback.parameters.maxIter
            if verbose
                displogger(lpvdd,'FB','|----------------------------- <strong>terminating : maxIter reached</strong> ----------------------------|\n');
            end
        end
    end
    prev_obj = obj.obj;
    
    
end % END FOR

if verbose
    t = toc;
    displogger(lpvdd,'FB','--- Finished (%s). Total run time: %f seconds.\n',char(datetime(now,'ConvertFrom','datenum')), t);
    
end

%% 
% Controller MAT OUT

K = lpvdd.FB.controller.K0mat;
lpvdd.Feedback.sol.obj = obj;
lpvdd.Feedback.sol.K0 = K;
end % END SOLVE
