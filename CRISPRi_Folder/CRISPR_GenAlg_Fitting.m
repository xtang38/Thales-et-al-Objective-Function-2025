clear all;
clc;
global par p

% Initial Values
par.Pytot = 0.5*10^(-9);  %%% 0.5nM constitutive GFP plasmid
par.dCas9tot=35*10^(-9);

%%% STAR Parameters: alpha_s, deg_s, deg_m, beta_s, alpha_m, KI, KE, alpha_gm

% par.alpha_cr, par.alpha_tr for transcription
Lalpha = 10^(-2);
Ualpha = 10;
stepRalpha = linspace(Lalpha,Ualpha,10);  
alpha_cr = randsample(stepRalpha,1);
alpha_tr = randsample(stepRalpha,1);

% par.deg_cr, par.deg_tr, par.deg_m for degradation
Ldeg = 10^(-5);
Udeg = 10^(-1);
stepDeg = linspace(Ldeg,Udeg,10); 
deg_cr = randsample(stepDeg,1);
deg_tr = randsample(stepDeg,1);
deg_m = randsample(stepDeg,1);

% par.omega for CRISPR binding
LoM = 10^3;
UoM = 10^7;
stepBeta = linspace(LoM,UoM,10); 
omega = randsample(stepBeta,1);

% par.gamma1 & par.gamma2 CRISPR formation rates
Lgam = 10^3;
Ugam = 10^7;
stepGam = linspace(Lgam,Ugam,10); 
gamma1 = randsample(stepGam,1);
gamma2 = randsample(stepGam,1);

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


p0 = [alpha_cr, alpha_tr, deg_cr, deg_tr, deg_m, omega, gamma1, gamma2, alpha_m, KI, KE, alpha_gm];

A = [];
b = [];
Aeq = [];
beq = [];
nlcon = [];

lb = [Lalpha Lalpha Ldeg Ldeg Ldeg LoM Lgam Lgam LaM LKi LKe LGm]; 
ub = [Ualpha Ualpha Udeg Udeg Udeg UoM Ugam Ugam UaM UKi UKe UGm]; 

Evalall = [];

options = optimoptions('ga','PlotFcn','gaplotbestf','FunctionTolerance',1e-8);

[p,fval] = ga(@CRISPR_GenAlg_Obj,length(p0),A,b,Aeq,beq,lb,ub,nlcon,options);

error = CRISPR_GenAlg_Obj(p);
% % savefile = strcat('CRISPR_GenAlg_Weighted_smooth2.mat');
% % save(savefile,"p0","p","error","fval")
disp(['Final SSE Objective: ' error])






