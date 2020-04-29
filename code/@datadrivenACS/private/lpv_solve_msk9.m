
function [K_out,obj_out] = lpv_solve_msk9(dd,W,n,verbose)
% DATA_DRIVEN SOLVE, USES MOSEK + MOSEK FUSION
import mosek.fusion.*;
Ts = dd.Feedback.controller.Ts;
z = tf('z',Ts);



M = Model('DATA-DRIVEN OPTIMIZATION');



slack = M.variable('slack',9,Domain.greaterThan(0));

if isempty(dd.FB.sol_obj)
    dd.FB.sol_obj.slack = 100;
end

if max(abs(dd.FB.sol_obj.slack))<1e-5
    % If last solution almost satisfies the constraints, fix them and
    % optimize over the objective
    M.constraint(slack,Domain.equalsTo(0));
    % l1 flag to optimize over constraints / objective
    if  dd.FB.internals.l1 == 1
        if verbose
        displogger(dd,'FB','|--------   FOUND CONTROLLER SATISFYING CONSTRAINTS - OPTIMIZING OVER OBJECTIVE   -------|\n');
        end
        dd.FB.internals.l1_hasChanged = 1;
        if dd.FB.internals.nIter == 1
            % if was possible to satisfy the constraints at the 1st
            % iteration, use the (previous) initial controller as starting
            % point. Avoid zero/pole near boundary.
           if verbose
            displogger(dd,'FB','|--- CONTROLLER SATISFYING CONSTRAINTS AT ITER 1 - REUSING INIT CONTROLLER FOR ITER 2 ---|\n');
           end
           formatK(dd);
        end
    end
    dd.FB.internals.l1 = 0;
end


%%
% Prepare initial controller
den = dd.FB.controller.K0mat.den;
num = dd.FB.controller.K0mat.num;

%{
Build the different powers of the scheduling parameters, ie 1, x, x^2, y,
y^2, xy, ... etc
%}
szy_prev = size(den,1);

Q = theta(dd);
schedParams = length(Q);


szx = 1+dd.Feedback.controller.order;
szy = 1+dd.Feedback.controller.order-order(1/dd.Feedback.controller.Fy);

if szy_prev > szy
    delta = szy_prev - szy;
    szx = szx + delta;
    szy = szy +delta;
else
    delta = 0;
    num = [zeros(-szy_prev + szy,length(Q));num];
    den = [zeros(-szy_prev + szy,length(Q));den];
    
end
sX_var = M.variable('X',[szx,schedParams]);
sY_var = M.variable('Y',[szy,schedParams]);


%%

%{

OBJECTIVE : 1st step find controller that satisfies the (soft) constraints,
then optimize over the objective using hard constraints.
%}
nCon = length(W);
nmod = dd.M.nmod;

% mmod -> multi-modoel
% if ~dd.FB.parameters.worstCase
%     nmod_ = nmod; % nmod_1 = 1 optimize for worst case
% else
    nmod_ = 1;
% end

gamma_2_mmod = M.variable('gamma2',[nCon,nmod_],Domain.greaterThan(0));


gamma_inf_mmod = M.variable('gamma_inf',[1,nmod_],Domain.greaterThan(0));
I = ones(nCon,1)/nmod_; 
I(n+1:end) = 0;

gamma_Inf_mmod = Expr.mul(Matrix.dense(ones(nCon,1)/nmod_),gamma_inf_mmod);


domega = diff([0,W])';
domega(n:end) = 0;
domega(1) = W(1);


norm_2 = M.variable('norm2',nCon,Domain.greaterThan(0));
M.constraint(Expr.sub(Expr.dot(norm_2,domega),slack.pick(-1+7)),Domain.lessThan(1));



domega = repmat(domega,1,nmod_);

obj_2 = M.variable('o_2',1);
M.constraint(Expr.sub(Expr.dot(domega,gamma_2_mmod),obj_2),Domain.lessThan(0));
%o = M.variable('o',Domain.equalsTo(1));
if dd.FB.internals.l1
    OBJ = Expr.sum(slack) ;
else
    OBJ = Expr.add([Expr.sum(Expr.mul(Matrix.dense(I),gamma_inf_mmod)), Expr.dot(Matrix.dense(domega/nmod_),gamma_2_mmod)]);
end

obj = M.variable('objective',1);
M.constraint(Expr.sub(OBJ,obj),Domain.lessThan(0)); % for obj.level();

M.objective('obj', ObjectiveSense.Minimize, OBJ);

%% ----------------

sa = size(den,1);
M.constraint(sY_var.index([sa,1]-1), Domain.greaterThan(0)); % Y(end) =1

for jj = 1 : delta
    for ii = 1 : (schedParams)
    M.constraint(sX_var.index([jj,ii]-1), Domain.equalsTo(0)); % X(jj) =0
    M.constraint(sY_var.index([jj,ii]-1), Domain.equalsTo(0)); % Y(jj) =0
    end
end
    


for ii = 2 : (schedParams)
    M.constraint(sY_var.index([sa,ii]-1), Domain.equalsTo(0)); % Y(end) =0
end

if dd.Feedback.parameters.lpv.fixedDC
    if length(Q) > 1
        for idx = 2 : length(Q)
            M.constraint(Expr.dot(ones([szx,1]),sX_var.slice([0,idx-1],[szx,idx])),Domain.equalsTo(0.0));
            M.constraint(Expr.dot(ones([szy,1]),sY_var.slice([0,idx-1],[szy,idx])),Domain.equalsTo(0.0));
        end
    end
end
%% ----------------

if ~isempty(dd.Feedback.parameters.c2)
        M.constraint(slack.pick(-1+9),Domain.lessThan(dd.Feedback.parameters.c2-1e-5));    
        % avoid change in sign
end
    
z_ = squeeze(freqresp(z,W));

Zy = (z_.^(0:(szy-1))); % [0,z,...,z^nx]
ZFy = Zy; 

if (~isempty(dd.FB.parameters.lpv.additionalPts) && dd.FB.parameters.lpv.enable)
    naddpts = size(dd.FB.parameters.lpv.additionalPts,1);
    for jj = 1 : naddpts
        Q = dd.FB.parameters.lpv.theta(dd.FB.parameters.lpv.additionalPts(jj,:));
        Q = Q(:);
        Y_c = (den*Q);
        Yc=Zy*(Y_c);
        x1_c = 2*real(conj(Yc).*ZFy)./abs(Yc);
        
        Y_n = Expr.mul(sY_var,Matrix.dense(Q));
        
        x1 = Expr.sub( Expr.mul(Matrix.dense(x1_c),Y_n), real(conj(Yc).*Yc)./abs(Yc));
        M.constraint(x1,Domain.greaterThan(dd.FB.parameters.c0));
    end
end

for mod =  1: nmod
    
    
    if nmod_ > 1
        gamma_2 = gamma_2_mmod.slice([0,mod-1],[nCon,nmod]);
        gamma_Inf = gamma_Inf_mmod.slice([0,mod-1],[nCon,nmod]);
    else
        gamma_Inf = gamma_Inf_mmod.slice([0,0],[nCon,1]);
        gamma_2 = gamma_2_mmod.slice([0,0],[nCon,1]);
    end
    
    % ----------------
    % Get correct scheluling parameter + compute K(Q) = X(Q)/Y(Q)
    Q = theta(dd,mod);
    if all(not(Q))
        error('At least one scheduling parameter must be non-zero')
    end
    X_n = Expr.mul(sX_var,Matrix.dense(Q));
    Y_n = Expr.mul(sY_var,Matrix.dense(Q));
    
    M.constraint(Y_n.index([szy,1]-1), Domain.equalsTo(1));
    % Y_var is monic
    
    
    X_c = (num*Q);
    Y_c = (den*Q);
    
    % ----------------
    if ~isempty(dd.Feedback.parameters.c2)
        M.constraint(Expr.add(...
        Expr.mul(sign(sum(X_c)),Expr.sum(X_n)), ... + 
        Expr.sub(slack.pick(-1+9),dd.Feedback.parameters.c2)),... 
        Domain.greaterThan(0)); % constraint at w=0, ie z = 1
    end
    
    XY_var = Expr.vstack(X_n,Y_n);
    

    
    %%
    
    z_ = squeeze(freqresp(z,W));
    Zy = (z_.^(0:szy-1))./(z_.^szy); % [0,z,...,z^nx]
    Zx = (z_.^(0:(szx-1)))./(z_.^szy); % [0,z,...,z^ny]
    
   
    G = dd.M.G(W,mod);% squeeze(freqresp(dd.G{mod},W))/dd.lambda;
    
    Fx = 1;%squeeze(freqresp(dd.controller.Fx,W));
    Fy = squeeze(freqresp(dd.Feedback.controller.Fy,W));
    
    Ycs=Zy*(Y_c);Xcs=Zx*(X_c);
    Yc = Ycs.*Fy;Xc = Xcs.*Fx;
    ZFy = Zy.*Fy; ZFx = Zx.*Fx;
    Pc = Yc + G.*Xc;
    
    
    x1_c = 2*real(conj(Zy).*Ycs)./abs(Ycs);
    x1 = Expr.sub(Expr.mul(Matrix.dense(x1_c),Y_n),real(Ycs.*conj(Ycs))./abs(Ycs));
    M.constraint(x1, Domain.greaterThan(dd.Feedback.parameters.c0));
    
    %  only do if controller satisifies constraints
    if ~dd.FB.internals.l1 %&& dd.Model.useForObjective(mod)
        
        Cp =  [G.*ZFx, ZFy].*conj(Pc); % P = Cp*[X;Y], Cp complex, [X;Y]  real
        
        
        lambdao1 =  (1./(min(abs(Pc),1))); % scaling factor
        x1_c = 2*real(Cp);
        x1 = Expr.sub( Expr.mul(Matrix.dense(x1_c.*lambdao1),XY_var), real(conj(Pc).*Pc.*lambdao1));
        
        % H_inf
        x2 = gamma_Inf;
        x3 = [];
        %
        if ~isempty(dd.FB.objective.oinfW1)
            oW1 = squeeze(dd.FB.objective.oinfW1(W,mod)).*sqrt(lambdao1);
            
            if any(abs(oW1)>1e8) && mod==1 && dd.flag.largeObj ~=1
                warning('Objective extremly large. Check frequencies');dd.flag.largeObj =1;
            end
            x3_d = oW1.*ZFy;
            x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
        end
        
        if ~isempty(dd.FB.objective.oinfW2)
            oW2 = squeeze(dd.FB.objective.oinfW2(W,mod)).*sqrt(lambdao1);
            if any(abs(oW2)>1e8) && mod==1 && dd.flag.largeObj ~=1
                warning('Objective extremly large. Check frequencies');dd.flag.largeObj =1;
            end
            x3_d = oW2.*G.*ZFx;
            if isempty(x3)
                x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            else
                x3 =  Expr.hstack( x3, Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            end
        end
        
        if ~isempty(dd.FB.objective.oinfW3)
            oW3 = squeeze(dd.FB.objective.oinfW3(W,mod)).*sqrt(lambdao1);
            if any(abs(oW3)>1e8) && mod==1 && dd.flag.largeObj ~=1
                warning('Objective extremly large. Check frequencies');dd.flag.largeObj =1;
            end
            x3_d = oW3.*ZFx;
            if isempty(x3)
                x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            else
                x3 =  Expr.hstack( x3, Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            end
        end
        
%         if ~isempty(dd.FB.objective.oinfW4)
%             oW4 = squeeze(dd.FB.objective.oinfW4(W,mod)).*sqrt(lambdao1);
%             x3_d = oW4.*G.*ZFy;
%             if isempty(x3)
%                 x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
%             else
%                 x3 =  Expr.hstack(x3, Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
%             end
%         end
        
        if ~isempty(x3)
            %M = RotatedConeToCone(M,(Expr.hstack(x1,x2,x3)));
            M.constraint((Expr.hstack(x1,x2,x3)), Domain.inRotatedQCone());
        end % END Hinf
        
        % H2
        
        x2 =  gamma_2;
        x3 = [];
        if ~isempty(dd.FB.objective.o2W1)
            oW1 = squeeze(dd.FB.objective.o2W1(W,mod)).*sqrt(lambdao1);
            x3_d = oW1.*ZFy;
            x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
        end
        
        if ~isempty(dd.FB.objective.o2W2)
            oW2 = squeeze(dd.FB.objective.o2W2(W,mod)).*sqrt(lambdao1);
            x3_d = oW2.*G.*ZFx;
            if isempty(x3)
                x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            else
                x3 =  Expr.hstack( x3, Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            end
        end
        
        if ~isempty(dd.FB.objective.o2W3)
            oW3 = squeeze(dd.FB.objective.o2W3(W,mod)).*sqrt(lambdao1);
            x3_d = oW3.*ZFx;
            if isempty(x3)
                x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            else
                x3 =  Expr.hstack( x3, Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
            end
        end
        
%         if ~isempty(dd.FB.objective.o2W4)
%             oW4 = squeeze(dd.FB.objective.o2W4(W,mod)).*sqrt(lambdao1);
%             x3_d = oW4.*G.*ZFy;
%             if isempty(x3)
%                 x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
%             else
%                 x3 =  Expr.hstack(x3, Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
%             end
%         end
        
        if ~isempty(x3)
            %M = RotatedConeToCone(M,(Expr.hstack(x1,x2,x3)));
            M.constraint((Expr.hstack(x1,x2,x3)), Domain.inRotatedQCone());
        end % END H2
        %end % END OBJECTIVE Loop-shaping not benched
    end
    % CONSTAINTS |W_n S_n|_inf< 1    W1, W2, W3, W4
    
    Cp =  [G.*ZFx, ZFy].*conj(Pc); % P = Cp*[X;Y], Cp complex, [X;Y]  real
    
    x1_c = 2*real(Cp);
    
    
    lambdaW1 =  (1./(abs(Pc).^2)); % scaling factor
    x1 = Expr.sub( Expr.mul(Matrix.dense(0.5*2*real(Cp).*lambdaW1),XY_var), 0.5*real(conj(Pc).*Pc).*lambdaW1 );
    
    
    if ~isempty(dd.FB.constraints.cinfW1) %&&  dd.Model.useForConstraint(mod)
        x2 = Expr.mul(Matrix.dense(ones(nCon,1)), Expr.add(slack.pick(-1+1),1));
        x3_d = (squeeze(dd.FB.constraints.cinfW1(W) )) .* ZFy .* sqrt(lambdaW1);
        x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
        M.constraint(Expr.hstack(x1,x2,x3), Domain.inRotatedQCone());
        %M = RotatedConeToCone(M,Expr.hstack(x1,x2,x3));
        
    end
    
    if (~isempty(dd.FB.constraints.cinfW2) || ~isempty(dd.FB.constraints.cinfW3)) %&&  dd.Model.useForConstraint(mod)
        if ~isempty(dd.FB.constraints.cinfW2)
            x3_d1 = abs(squeeze(dd.FB.constraints.cinfW2(W)).*G);
        else
            x3_d1 = zeros(nCon,1);
        end
        
        if ~isempty(dd.FB.constraints.cinfW3)
            x3_d2 = abs(squeeze(dd.FB.constraints.cinfW3(W)));
        else
            x3_d2 = zeros(nCon,1);
        end
        
        x3_d2 = 1./max(x3_d1,x3_d2).^2;
        x3_d = ZFx .* sqrt(lambdaW1);
        x2 = Expr.mul(Matrix.dense(x3_d2), Expr.add(slack.pick(-1+3),1));
        
        x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
        M.constraint((Expr.hstack(x1,x2,x3)), Domain.inRotatedQCone([nCon,4]));
        %M = RotatedConeToCone(M,(Expr.hstack(x1_scaled,x2_scaled,x3)));
        
    end
    
%     if ~isempty(dd.FB.constraints.cinfW4) %&&  dd.Model.useForConstraint(mod)
%         W4 = squeeze(dd.FB.constraints.cinfW4(W));
%         x1 = 2.*real(conj(Yc).*ZFy);
%         z1 = Expr.sub( Expr.mul( Matrix.dense(x1.*lambdaW1),Y_n), Matrix.dense(real(conj(Yc).*Yc).*lambdaW1));
%         z2 = Expr.mul(Matrix.dense(0.5*ones(nCon,1)), Expr.add(slack.pick(-1+4),1));
%         x3_d = W4.*ZFx.* sqrt(lambdaW1);
%         z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
%         M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone([nCon,4]));
% 
%     end
        
    % CONSTAINTS |W_n S_n|_2< 1    W1, W2, W3, W4
    % Philippe Schuchert - Oct. 2019. 
    % norm 2 constrants could easly be
    % implemented, but there is no need for them. When a 2 norm is usually
    % used, it is in the objective. not benched yet.
   
%     if ~isempty(dd.FB.constraints.c2W1)
%          lambdao1 = abs(squeeze(dd.FB.constraints.c2W1(W)));
%          x1 = 2.*real(conj(Yc).*ZFy);
%         W3 = squeeze(dd.FB.constraints.c2W3(W)).*sqrt(lambdao1);
%         z1 = Expr.sub( Expr.mul( Matrix.dense(x1.*lambdaW1),Y_n), Matrix.dense(real(conj(Yc).*Yc).*lambdaW1));
%         
%         x2 =  norm_2;
%         x3_d = W3.*ZFy.* sqrt(lambdaW1);
%         x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),Y_n),Expr.mul(Matrix.dense(imag(x3_d)),Y_n));
%         M.constraint((Expr.hstack(z1,x2,x3)), Domain.inRotatedQCone());
%     end
%     
%     if ~isempty(dd.FB.constraints.c2W2)
%         lambdao1 = abs(squeeze(dd.FB.constraints.c2W2(W)));
%          x1 = 2.*real(conj(Yc).*ZFy);
%         W3 = sqrt(squeeze(abs(dd.FB.constraints.c2W3(W)))).*sqrt(lambdao1);
%         z1 = Expr.sub( Expr.mul( Matrix.dense(x1.*lambdaW1),Y_n), Matrix.dense(real(conj(Yc).*Yc).*lambdaW1));
%         
%         x2 =  norm_2;
%         x3_d = W3.*ZFx.* sqrt(lambdaW1);
%         x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
%         M.constraint((Expr.hstack(z1,x2,x3)), Domain.inRotatedQCone());
%     end
%     
%     if ~isempty(dd.FB.constraints.c2W3)
%         
%         W3 = (abs(squeeze(dd.FB.constraints.c2W3(W))));
%         
%         lambdao1 =  1./abs(max(0.001,abs(W3.*Pc.*conj(Pc)))); % scaling factor
%         x1_c = 2*real(Cp);
%         z1 = Expr.sub( Expr.mul(Matrix.dense(x1_c.*lambdao1),XY_var), real(conj(Pc).*Pc).*lambdao1);
%         x2 =  norm_2;      
% 
%         x3_d = W3.*ZFx.* sqrt(lambdaW1);
%         x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
%         M.constraint((Expr.hstack(z1,x2,x3)), Domain.inRotatedQCone());
%     end
    
%     if ~isempty(dd.FB.constraints.c2W4)
%         lambdao1 = abs(squeeze(dd.FB.constraints.c2W4(W)));
%         x1 = 2.*real(conj(Yc).*ZFy);
%         W3 = squeeze(dd.FB.constraints.c2W3(W)).*sqrt(lambdao1);
%         z1 = Expr.sub( Expr.mul( Matrix.dense(x1.*lambdaW1),Y_n), Matrix.dense(real(conj(Yc).*Yc).*lambdaW1));
%         
%         x2 =  Expr.dot(Matrix.dense(sqrt(2)),norm_2);
%         x3_d = W3.*ZFx.* sqrt(lambdaW1);
%         x3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_n),Expr.mul(Matrix.dense(imag(x3_d)),X_n));
%         M.constraint((Expr.hstack(z1,x2,x3)), Domain.inRotatedQCone());
%     end
    
    
    % ----- 
    M.constraint(Expr.mul(Matrix.dense(x1_c./abs(Pc)),XY_var),Domain.greaterThan(dd.Feedback.parameters.c1)); 
end
M.setSolverParam('intpntCoTolRelGap', 1e-6)
M.setSolverParam('intpntCoTolMuRed', 1e-6)


M.solve();
M.acceptedSolutionStatus(AccSolutionStatus.Anything);

kx = sX_var.level();
ky = sY_var.level();

ZFx = reshape(kx,[schedParams, szx])';
ZFy = reshape(ky,[schedParams, szy])';



ZFx(1:delta,:) = [];
ZFy(1:delta,:) = [];

K_out.num = ZFx;
K_out.den = ZFy;
K_out.Ts = Ts;


obj_out.H2        = sqrt(obj_2.level());
obj_out.Hinf      = sqrt(max(gamma_inf_mmod.level())/(ZFy(end,1)))*sqrt(2);
obj_out.slack     = max(abs(slack.level()));
obj_out.primal    = char(M.getPrimalSolutionStatus);
obj_out.dual      = char(M.getDualSolutionStatus);
obj_out.primalVal = M.primalObjValue();
obj_out.dualVal   = M.dualObjValue();
obj_out.obj       = max(obj.level(),0);

M.dispose();

end % END DATA_DRIVEN_SOLVE

