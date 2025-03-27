clear all;
clc;
global par p iter


% Initial Values
par.PYtot = 0.5*10^(-9);  %%% 0.5nM STAR targeting GFP plasmid

%%% STAR Parameters: alpha_s, deg_s, deg_m, beta_s, alpha_m, KI, KE, alpha_gm


% par.alpha_s for transcription
Lalpha = 10^(-2);
Ualpha = 10;
stepRalpha = linspace(Lalpha,Ualpha,10);  
alpha_s = randsample(stepRalpha,1);

% par.deg_s & par.deg_m for degradation
Ldeg = 10^(-5);
Udeg = 10^(-1);
stepDeg = linspace(Ldeg,Udeg,10); 
deg_s = randsample(stepDeg,1);
deg_m = randsample(stepDeg,1);

% par.beta for STAR binding
Lbd = 10^3;
Ubd = 10^7;
stepBeta = linspace(Lbd,Ubd,10); 
beta = randsample(stepBeta,1);

% par.alpha_m for GFP mRNA transcription rate
LaM = 10^(-2);
UaM = 100;
stepRaM = linspace(LaM,UaM,10);
alpha_m = randsample(stepRaM,1);

% par.KI for initiation rate
LKi = 10^(-4);
UKi = 10^(-2);
stepKi = linspace(LKi,UKi,10); 
KI = randsample(stepKi,1);

% par.KE for elongation rate 
LKe = 10^(-4);
UKe = 10^(-2);
stepKe=linspace(LKe,UKe,10);
KE = randsample(stepKe,1);

% par.alpha_gm for maturation rate
LGm = 10^(-3);
UGm = 10^(-1);
stepGm = linspace(LGm,UGm,10); 
alpha_gm = randsample(stepGm,1);


for iter = 1:100  
    p0 = [alpha_s, deg_s, deg_m, beta, alpha_m, KI, KE, alpha_gm];
    
    A = [];
    b = [];
    Aeq = [];
    beq = [];
    nlcon = [];

    lb = [Lalpha Ldeg Ldeg Lbd LaM LKi LKe LGm]; 
    ub = [Ualpha Udeg Udeg Ubd UaM UKi UKe UGm]; 
    
    Evalall = [];
    
    options = optimoptions('ga','PlotFcn','gaplotbestf','FunctionTolerance',1e-8);
    
    [p,fval] = ga(@STAR_GenAlg_Obj,length(p0),A,b,Aeq,beq,lb,ub,nlcon,options);
    
    

    error = STAR_GenAlg_Obj(p);
    % % savefile = strcat('STAR_syndata_smooth_weighted',int2str(iter),'.mat');
    % % save(savefile,"p0","p","error","fval")
    
end





