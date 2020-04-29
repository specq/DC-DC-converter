function [K_out,obj_out] = lpv_solve_ff_msk9(dd,W)
% DATA_DRIVEN SOLVE, USES MOSEK + MOSEK FUSION


Ts = dd.FF.controller.K0mat.Ts;
z = tf('z',Ts);


import mosek.fusion.*;
M = Model('DATA-DRIVEN OPTIMIZATION');

%%
try
    KLPV = computeK(dd);
catch
    % FB not computed
    KLPV = tf(0);
end
% Prepare initial controller
den = dd.FF.controller.K0mat.den;
num = dd.FF.controller.K0mat.num;

%{
Build the different powers of the scheduling parameters, ie 1, x, x^2, y,
y^2, xy, ... etc
%}
theta_ = theta(dd,1,'FF');
schedParams = length(theta_);
sX_var = M.variable('X',[1+dd.FF.controller.order,schedParams]);
sY_var = M.variable('Y',[1+dd.FF.controller.order,schedParams]);

%%


nCon = length(W);
nmod = dd.M.nmod;
nmod_ = dd.M.nmod;
% mmod -> multi-modoel
gamma_2_mmod = M.variable('gamma2',[nCon,nmod_],Domain.greaterThan(0));
gamma_inf_mmod = M.variable('gamma_inf',Domain.greaterThan(0));
gamma_Inf_mmod = Expr.mul(Matrix.dense(ones(nCon,nmod)),gamma_inf_mmod);

%norm_2 = M.variable('norm2',nCon,Domain.greaterThan(0));
%M.constraint(Expr.dot(norm_2,domega),Domain.lessThan(1));
domega = diff([0,W])';


%


slack = M.variable('slack',8,Domain.greaterThan(0));
if isempty(dd.FF.sol_obj)
    dd.FB.sol_obj.slack = 100;
end

dd.internals.l1 = 0;

% Objective

domega = repmat(domega,1,nmod_);

obj_2 = M.variable('o_2',1,Domain.greaterThan(0));
M.constraint(Expr.sub(Expr.dot((domega),gamma_2_mmod),obj_2),Domain.equalsTo(0));

OBJ = Expr.add(Expr.sum(gamma_inf_mmod), obj_2);


obj = M.variable('objective',1);
M.constraint(Expr.sub(OBJ,obj),Domain.equalsTo(0)); % for obj.level();

M.objective('obj', ObjectiveSense.Minimize, OBJ);

%% ----------------

sa = size(den,1);
M.constraint(sY_var.index([sa,1]-1), Domain.greaterThan(0)); % Y(end) =1
for ii = 2 : (schedParams)
    M.constraint(sY_var.index([sa,ii]-1), Domain.equalsTo(0)); % Y(end) =1
end

Q = theta(dd,1,'FF');
if dd.FF.parameters.fixedDC
    if length(Q) > 1
        for idx = 2 : length(Q)
            M.constraint(Expr.dot(ones([1+dd.FF.controller.order,1]),sX_var.slice([0,idx-1],[dd.FF.controller.order+1,idx])),Domain.equalsTo(0.0));
            M.constraint(Expr.dot(ones([1+dd.FF.controller.order,1]),sY_var.slice([0,idx-1],[dd.FF.controller.order+1,idx])),Domain.equalsTo(0.0));
            
        end
    end
end


% 
% if dd.FF.parameters.fixedHF
%     if length(Q) > 1
%         for idx = 2 : length(Q)
%             O = ones([1+dd.FF.controller.order,1]);
%             O(2:2:end) = -1;
%             M.constraint(Expr.dot(O,sX_var.slice([0,idx-1],[dd.FF.controller.order+1,idx])),Domain.equalsTo(0.0));
%             M.constraint(Expr.dot(O,sY_var.slice([0,idx-1],[dd.FF.controller.order+1,idx])),Domain.equalsTo(0.0));
%         end
%     end
% end
z_ = squeeze(freqresp(z,W));

Zy = (z_.^(0:(dd.FF.controller.order))); % [0,z,...,z^nx]
Zx = (z_.^(0:(dd.FF.controller.order))); % [0,z,...,z^ny]

Ky = Zy; Kx = Zx;


if (~isempty(dd.FF.parameters.lpv.additionalPts) && dd.FF.parameters.lpv.enable)
    naddpts = size(dd.FF.parameters.lpv.additionalPts,1);
    for jj = 1 : naddpts
        Q = dd.FF.parameters.lpv.theta(dd.FF.parameters.lpv.additionalPts(jj,:));
        Q = Q(:);
        Y_c = (den*Q);
        Yc=Zy*(Y_c);
        x1_c = 2*real(conj(Yc).*Ky)./abs(Yc);
        
        Y_var = Expr.mul(sY_var,Matrix.dense(Q));
        
        x1 = Expr.sub( Expr.mul(Matrix.dense(x1_c),Y_var), real(conj(Yc).*Yc)./abs(Yc));
        M.constraint(x1,Domain.greaterThan(dd.FF.parameters.c0));
    end
end
for mod =  1: nmod
    if length(KLPV)>1
        KLPV_loc = squeeze(freqresp(KLPV(:,:,mod),W));
    else
        KLPV_loc = squeeze(freqresp(KLPV,W));
    end
    if nmod_ > 1
        gamma_2 = gamma_2_mmod.slice([0,mod-1],[nCon,mod]);
    else
        gamma_2 = gamma_2_mmod.slice([0,0],[nCon,1]);
    end
    gamma_Inf = gamma_Inf_mmod.slice([1,nmod]-1,[nCon+1,nmod+1]-1);
    
    % ----------------
    % Get correct scheluling parameter + compute K(Q) = X(Q)/Y(Q)
    Q = theta(dd,mod,'FF');
    
    if all(not(Q))
        error('At least one scheduling parameter must be non-zero')
    end
    X_var = Expr.mul(sX_var,Matrix.dense(Q));
    Y_var = Expr.mul(sY_var,Matrix.dense(Q));
    XY_var = Expr.vstack(X_var,Y_var);
    
    
    
    
    M.constraint(Y_var.index([1+dd.FF.controller.order,1]-1), Domain.equalsTo(1));
    
    X_c = (num*Q);
    Y_c = (den*Q);
    
    Yc=Zy*(Y_c);
    Xc=Zx*(X_c);
    
    
    G = dd.M.G(W,mod);% 
    
    
    if 1 % ~dd.internals.l1
        x1_c = 2*real(conj(Yc).*Ky)./abs(Yc);
        x1 = Expr.sub( Expr.mul(Matrix.dense(x1_c),Y_var), real(conj(Yc).*Yc)./abs(Yc));
        M.constraint(x1,Domain.greaterThan(dd.FF.parameters.c0));
        
        
        
        if ~isempty(dd.FF.objective.o2W1) && dd.Model.useForObjective(mod)
            oW1 = (squeeze(dd.FF.objective.o2W1(W)));
            alpha = (1+G.*KLPV_loc);
            
            lambda =  1./abs(alpha);
            x3_d = [-G.*Kx,Ky]./alpha.*sqrt(lambda).*oW1;
            z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),XY_var),Expr.mul(Matrix.dense(imag(x3_d)),XY_var));
            
            z2 = gamma_2;
            
            
            x1 = 2.*real(conj(Yc).*Ky);
            z1 = Expr.sub( Expr.mul(Matrix.dense( x1.*lambda),Y_var), real(conj(Yc).*Yc).*lambda);
            
            M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone());
        end
        
        if ~isempty(dd.FF.objective.oinfW1) && dd.Model.useForObjective(mod)
            oW1 = (squeeze(dd.FF.objective.oinfW1(W)));
            alpha = (1+G.*KLPV_loc);
            
            lambda = 1./abs(alpha);
            x3_d = [-G.*Kx,Ky]./alpha.*(lambda).*oW1;
            z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),XY_var),Expr.mul(Matrix.dense(imag(x3_d)),XY_var));
            
            z2 = Expr.mulElm(Matrix.dense(0.5*ones(size(G)).*lambda),gamma_Inf);
            
            
            x1 = 2.*real(conj(Yc).*Ky);
            z1 = Expr.sub( Expr.mul( Matrix.dense(x1.*lambda),Y_var), real(conj(Yc).*Yc).*lambda);
            
            % M.constraint(z1,Domain.greaterThan(dd.c0));
            M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone());
        end
        
        if ~isempty(dd.FF.objective.oinfW3) && dd.Model.useForObjective(nmod)
            alpha = 1;%(1+G.*KLPV_loc);
            
            lambda = 1./abs(Yc);
            x1 = 2.*real(conj(Yc).*Ky);
            z1 = Expr.sub( Expr.mul( x1.*lambda,Y_var), real(conj(Yc).*Yc.*lambda));
            ffcinfW3 = squeeze(dd.FF.objective.oinfW3(W));
            x3_d = [-Kx,KLPV_loc.*Ky].*(ffcinfW3).*sqrt(lambda)./alpha;
            z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),XY_var),Expr.mul(Matrix.dense(imag(x3_d)),XY_var));
            
            z2 =gamma_Inf;
            
            M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone());
        end
        
        
        if ~isempty(dd.FF.constraints.cinfW1) && dd.Model.useForConstraint(mod)
            
            cinfW1 = squeeze(dd.FF.constraints.cinfW1(W));
            
            alpha = (1+G.*KLPV_loc);
            
            
            x3_d = [-G.*Kx,Ky].*abs(cinfW1)./alpha;
            z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),XY_var),Expr.mul(Matrix.dense(imag(x3_d)),XY_var));
            
            z2 = (Expr.constTerm(0.5*ones(nCon,1)));
            
            
            x1 = 2.*real(conj(Yc).*Ky);
            z1 = Expr.sub( Expr.mul( Matrix.dense(x1),Y_var), real(conj(Yc).*Yc));
            
            M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone());
        end
        
        if ~isempty(dd.FF.constraints.cinfW2) && dd.Model.useForConstraint(mod)
            
            cinfW2= squeeze(dd.FF.constraints.cinfW2(W));
            alpha = (1+G.*KLPV_loc)./G;
            
             x3_d = [Kx,KLPV_loc.*Ky].*abs(cinfW2)./alpha;
            z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),XY_var),Expr.mul(Matrix.dense(imag(x3_d)),XY_var));
            
            z2 = (Expr.constTerm(0.5*ones(nCon,1)));
            
            
            x1 = 2.*real(conj(Yc).*Ky);
            z1 = Expr.sub( Expr.mul( Matrix.dense(x1),Y_var), real(conj(Yc).*Yc));
            
            M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone());
        end
        
        if ~isempty(dd.FF.constraints.cinfW3) && dd.Model.useForConstraint(mod)
            
            cinfW3= squeeze(dd.FF.constraints.cinfW3(W));
            alpha = (1+G.*KLPV_loc);
            
            x3_d = [Kx,KLPV_loc.*Ky].*abs(cinfW3)./alpha;
            z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),XY_var),Expr.mul(Matrix.dense(imag(x3_d)),XY_var));
            
            z2 = (Expr.constTerm(0.5*ones(nCon,1)));
            
            
            x1 = 2.*real(conj(Yc).*Ky);
            z1 = Expr.sub( Expr.mul( Matrix.dense(x1),Y_var), real(conj(Yc).*Yc));
            
            M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone());
        end
        
        
        if ~isempty(dd.FF.constraints.cinfW4) && dd.Model.useForConstraint(mod)
            
            cinfW4= squeeze(dd.FF.constraints.cinfW4(W));
            
            
            x3_d = Kx.*cinfW4;
            z3 =  Expr.hstack( Expr.mul(Matrix.dense(real(x3_d)),X_var),Expr.mul(Matrix.dense(imag(x3_d)),X_var));
            
            z2 = (Expr.constTerm(0.5*ones(nCon,1)));
            
            
            x1 = 2.*real(conj(Yc).*Ky);
            z1 = Expr.sub( Expr.mul( Matrix.dense(x1),Y_var), real(conj(Yc).*Yc));
            
            M.constraint((Expr.hstack(z1,z2,z3)), Domain.inRotatedQCone());
        end
    end
end
M.setSolverParam('intpntCoTolRelGap', 1e-6)
M.setSolverParam('intpntCoTolMuRed', 1e-6)

M.solve();
M.acceptedSolutionStatus(AccSolutionStatus.Anything);

kx = sX_var.level();
ky = sY_var.level();

Kx = reshape(kx,[schedParams, 1+dd.FF.controller.order])';
Ky = reshape(ky,[schedParams, 1+dd.FF.controller.order])';




K_out.num = Kx;
K_out.den = Ky;
K_out.Ts  = Ts;




obj_out.H2        = sqrt(obj_2.level());
obj_out.Hinf      = sqrt(max(gamma_inf_mmod.level())/(Ky(end,1)));
obj_out.slack     = max(abs(slack.level()));
obj_out.primal    = char(M.getPrimalSolutionStatus);
obj_out.dual      = char(M.getDualSolutionStatus);
obj_out.primalVal = M.primalObjValue();
obj_out.dualVal   = M.dualObjValue();
obj_out.obj       = max(obj.level(),0);

M.dispose();

end % END DATA_DRIVEN_SOLVE

