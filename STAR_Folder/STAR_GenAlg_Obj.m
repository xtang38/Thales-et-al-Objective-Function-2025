function error = STAR_GenAlg_Obj(p) 

global par iter

load Simulated_Data_Jan22_2025.mat;


%%% getting the stablized variance data
Var_weight = readmatrix("STAR_Smoothed_Variance.xlsx");
for ij = 1:length(Var_weight(:,1))
    maximum = max(Var_weight(ij,:));
    Var_weight(ij,:) = Var_weight(ij,:)./maximum;
end


tspan=0:300:14400;
x0=[0 0 0 0 0];

% Relative error & Abs tolerance for ODE Solver
options = odeset('RelTol',1e-10,'AbsTol',1e-10);    

par.Ps=8*10^(-9); %%% STAR 8nM

[t,x]=ode23s(@(t,x) STAR_GenAlg_Model(t,x,p),tspan,x0,options);

Data = Simu_Traje((iter-1)*3+1:iter*3,:);
Mean_Data = mean(Data);

GFP = x(:,5)*10^6;
Var8nM = var(Data)./max(var(Data));


error = 0;
%%%% Weighted with stablized variance

for i = 1:length(GFP)-1  %%%% first time point is omitted as it starts from 0.
    error = error + (GFP(i+1) - Mean_Data(i+1))^2/Var_weight(iter,i+1); 
end

%%%% Weighted with raw variance

% for i = 1:length(Var8nM)-1  %%%% first time point is omitted as it starts from 0.
%     error = error + (GFP(i+1) - Mean_Data(i+1))^2/Var8nM(i+1);
% end

%%%% unweighted
%%%% error = sum((GFP - Mean_Data').^2);  

clear GFP Data Mean_Data


end