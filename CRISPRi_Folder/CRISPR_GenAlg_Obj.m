function error = CRISPR_GenAlg_Obj(p) 

global par 

% Average values of 3 CRISPR experiments
M09=csvread('NewData3009.csv',15,0);


Var_raw = readmatrix("CRIPSRi_Smoothed_Variance.csv");
Var_weight = Var_raw(2:end,2);
Var_weight = Var_weight./max(Var_weight);


% Simulation timespan
tspan=0:300:14400;

par.Pytot = 0.5*10^(-9);
par.dCas9tot =35*10^(-9);
x0=[0 0 0 par.dCas9tot 0 0 0 0 0 0];

% Relative error & Abs tolerance for ODE Solver
options = odeset('RelTol',1e-10,'AbsTol',1e-10);    
   

%%% 0.25 nM CRISPR
par.Pcr=0.25*10^(-9);
par.Ptr=0.25*10^(-9);
        
         
[t,x]=ode23s(@(t,x) CRISPR_GenAlg_Model(t,x,p),tspan,x0,options); 

Data = M09(:,20:28)';
Mean_Data = mean(Data);

GFP = x(:,10)*10^6;
Var8nM = var(Data)./max(var(Data));


error = 0;

%%% WtS
for i = 1:length(GFP)-1  %%%% first time point is omitted as it starts from 0.
    error = error + (GFP(i+1) - Mean_Data(i+1))^2/Var_weight(i+1);  
end


%%% Wt
% % for i = 1:length(Var8nM)-1  %%%% first time point is omitted as it starts from 0.
% %     error = error + (GFP(i+1) - Mean_Data(i+1))^2/Var8nM(i+1);
% % end

 %%%% UNweighted
% % for i = 1:length(Var8nM)-1  %%%% first time point is omitted as it starts from 0.
% %     error = error + (GFP(i+1) - Mean_Data(i+1))^2;
% % end


clear GFP Data Mean_Data


end