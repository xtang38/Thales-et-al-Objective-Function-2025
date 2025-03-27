%%% simulated data generation with TRUE parameter fitted from STAR par.Ps = 8*10^(-9).
%%%

clear all;
clc;
close all
global par p

%%% Load TRUE kinetic parameters  %%%
load STAR_Nominal_Parameter.mat;
p=(p);

% Initial Values
par.PYtot = 0.5*10^(-9);  %%% 0.5nM STAR targeting GFP plasmid

par.Ps = 8*10^(-9); %%% STAR 8nM
x0 = [0 0 0 0 0];
tspan=0:300:14400; %%% seconds
options = odeset('RelTol',1e-10,'AbsTol',1e-10);
[t,x] = ode23s(@(t,x) STAR_GenAlg_Model(t,x,p),tspan,x0, options);
x = x.*(10^6);
Simu_t = t./60;
TRUE_Simu = x(:,5);

figure
plot(Simu_t,TRUE_Simu,'LineWidth',2)
title('STAR=4nM ')
xlabel('Time (min)')
ylabel('EGFP Conc. (\muM)')
xlim([0 250])
% ylim([0 0.1])
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')




%%%% Gaussian distribution data generation

L = length(tspan);

mu = linspace(0,max(TRUE_Simu)*0.2,L);
sigma = mu*0.2;

Trial = 300;
Normal_noise = zeros(Trial,L);

for ij = 1:Trial
    for i = 1:L
        Normal_noise(ij,i) = random('Normal',0,sigma(i));
    end
end

Simu_Traje = Normal_noise + TRUE_Simu';

figure
for ij = 1:Trial
    plot(Simu_t,Simu_Traje(ij,:),'k','LineWidth',1)
    hold on 
end
plot(Simu_t,TRUE_Simu,'-o','LineWidth',2)
hold off
title('STAR=8nM Simu_Individual')
xlabel('Time (min)')
ylabel('EGFP Conc. (\muM)')
xlim([0 250])
% ylim([0 0.06])
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')

Simu_Ave = zeros(Trial/3,L);
for j = 1:Trial/3
    Simu_Ave(j,:) = mean(Simu_Traje((j-1)*3+1:j*3,:));
end

figure
for ij = 1:Trial/3
    plot(Simu_t,Simu_Ave(ij,:),'k','LineWidth',1)
    hold on 
end

plot(Simu_t,TRUE_Simu,'-o','LineWidth',2)
hold off
title('STAR=8nM Simu Ave')
xlabel('Time (min)')
ylabel('EGFP Conc. (\muM)')
xlim([0 250])
% ylim([0 0.06])
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')

save('Simulated_Data_Jan22_2025.mat','Simu_Ave','Simu_Traje','Simu_t')


