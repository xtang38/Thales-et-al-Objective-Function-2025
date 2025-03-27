clear all;
clc;
close all

global par p

%%% Load Results %%%
    
% Initial Values
par.Pytot = 0.5*10^(-9);  %%% 0.5nM STAR targeting GFP plasmid
par.dCas9tot = 35*10^(-9);

Simu=cell(3,1);

for ij = 1:3
    if ij == 1
        load CRISPR_GenAlg_Weighted.mat;
       
        for i = 1:3
            if i==1  
                %%% 0.1 CRISPR
                par.Pcr=0.1*10^(-9);
                par.Ptr=0.1*10^(-9);
                
            elseif i==2
                %%% 0.25 nM CRISPR
                par.Pcr=0.25*10^(-9);
                par.Ptr=0.25*10^(-9);
                
            else
                %%% 0.5 nM CRISPR
                par.Pcr=0.5*10^(-9);
                par.Ptr=0.5*10^(-9);
            end
        
            x0 = [0 0 0 par.dCas9tot 0 0 0 0 0 0];
            tspan=0:300:14400; %%% seconds
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
            [t,x] = ode23s(@(t,x) CRISPR_GenAlg_Model(t,x,p),tspan,x0, options);
            x = x.*(10^6);
            Simu_t = t./60;
            Simu{ij,1}(:,i)=x(:,10);
          end
             clear p x

    elseif ij == 2
        load CRISPR_GenAlg_Weighted_smooth.mat;
        for i = 1:3
            if i==1  
                %%% 0.1 CRISPR
                par.Pcr=0.1*10^(-9);
                par.Ptr=0.1*10^(-9);
                
            elseif i==2
                %%% 0.25 nM CRISPR
                par.Pcr=0.25*10^(-9);
                par.Ptr=0.25*10^(-9);
                
            else
                %%% 0.5 nM CRISPR
                par.Pcr=0.5*10^(-9);
                par.Ptr=0.5*10^(-9);
            end
        
            x0 = [0 0 0 par.dCas9tot 0 0 0 0 0 0];
            tspan=0:300:14400; %%% seconds
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
            [t,x] = ode23s(@(t,x) CRISPR_GenAlg_Model(t,x,p),tspan,x0, options);
            x = x.*(10^6);
            Simu_t = t./60;
            Simu{ij,1}(:,i)=x(:,10);
          end
             clear p x

    else
        load CRISPR_GenAlg_Unweighted.mat;
        
        for i = 1:3
            if i==1  
                %%% 0.1 CRISPR
                par.Pcr=0.1*10^(-9);
                par.Ptr=0.1*10^(-9);
                
            elseif i==2
                %%% 0.25 nM CRISPR
                par.Pcr=0.25*10^(-9);
                par.Ptr=0.25*10^(-9);
                
            else
                %%% 0.5 nM CRISPR
                par.Pcr=0.5*10^(-9);
                par.Ptr=0.5*10^(-9);
            end
        
            x0 = [0 0 0 par.dCas9tot 0 0 0 0 0 0];
            tspan=0:300:14400; %%% seconds
            options = odeset('RelTol',1e-10,'AbsTol',1e-10);
            [t,x] = ode23s(@(t,x) CRISPR_GenAlg_Model(t,x,p),tspan,x0, options);
            x = x.*(10^6);
            Simu_t = t./60;
            Simu{ij,1}(:,i)=x(:,10);
          end
             clear p x

    end
    

end

M11=csvread('NewData3009.csv',15,0);
Time = M11(:,1); %%% mins

Data01 = M11(:,11:19);
Data025 = M11(:,20:28);  %%% fitted to 0.25 nM case
Data05 = M11(:,29:37);

Matrix01 = zeros(3,49);
Matrix025 = zeros(3,49);
Matrix05 = zeros(3,49);


for i = 1:49

    Matrix01(1,i) = min(Data01(i,:));  %%% 1 is for 4nm Weighted
    Matrix01(2,i) = mean(Data01(i,:)); %%% 2 is for 4nm Weighted stable
    Matrix01(3,i) = max(Data01(i,:)); %%% 3 is for 4nm unWeighted

    Matrix025(1,i) = min(Data025(i,:));  %%% 1 is for 4nm Weighted
    Matrix025(2,i) = mean(Data025(i,:));
    Matrix025(3,i) = max(Data025(i,:));

    Matrix05(1,i) = min(Data05(i,:));  %%% 1 is for 4nm Weighted
    Matrix05(2,i) = mean(Data05(i,:));
    Matrix05(3,i) = max(Data05(i,:));

end


%%%% green is unweighted; red is weighted
cooo = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
col =["[0 0.4470 0.7410]";"[0.8500 0.3250 0.0980]";"[0.9290 0.6940 0.1250]"];

colors = [0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];

Hcolors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0 0.4470 0.7410];


sc = 2;

figure
subplot(3,1,1) %%%% fitted comparison

plot(Time,Matrix025(2,:),'k','LineWidth',sc)
hold on
tx = [Time',fliplr(Time')];
py = [Matrix025(1,:),fliplr(Matrix025(3,:))];
fill(tx,py,'k','FaceAlpha',0.1,'EdgeColor','none');
hold on

plot(Time,Simu{3,1}(:,2),'Color',colors(1,:),'LineWidth',sc)  %%% unweighted
hold on
plot(Time,Simu{1,1}(:,2),'Color',colors(2,:),'LineWidth',sc)  %%% weighted
hold on
plot(Time,Simu{2,1}(:,2),'Color',colors(3,:),'LineWidth',sc)  %%% WtS
hold off
%xlabel('Time (mins)')
xlim([0 240])
% ylabel('EGFP Conc. (\muM) - 0.25nM')
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')

subplot(3,1,2)  %%%% prediction for 0.1nM

plot(Time,Matrix01(2,:),'k','LineWidth',sc)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [Matrix01(1,:),fliplr(Matrix01(3,:))];
fill(tx,py,'k','FaceAlpha',0.1,'EdgeColor','none');
hold on
plot(Time,Simu{3,1}(:,1),'Color',colors(1,:),'LineWidth',sc)  %%% unweighted
hold on
plot(Time,Simu{1,1}(:,1),'Color',colors(2,:),'LineWidth',sc)  %%% weighted
hold on
plot(Time,Simu{2,1}(:,1),'Color',colors(3,:),'LineWidth',sc)  %%% WtS
hold off
%xlabel('Time (mins)')
ylabel('EGFP Conc. (\muM)')
xlim([0 240])
set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')



subplot(3,1,3)  %%%% prediction for 0.1nM

plot(Time,Matrix05(2,:),'k','LineWidth',sc)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [Matrix05(1,:),fliplr(Matrix05(3,:))];
fill(tx,py,'k','FaceAlpha',0.1,'EdgeColor','none');
hold on
plot(Time,Simu{3,1}(:,3),'Color',colors(1,:),'LineWidth',sc)  %%% unweighted
hold on
plot(Time,Simu{1,1}(:,3),'Color',colors(2,:),'LineWidth',sc)  %%% weighted
hold on
plot(Time,Simu{2,1}(:,3),'Color',colors(3,:),'LineWidth',sc)  %%% WtS
hold off

xlim([0 240])
% % ylim([0 0.115])
% ylabel('EGFP Conc. (\muM) - 0.5nM')

% legend('Ave. Exp.','Distri. Exp.',...
    % 'Var Weighted','VarSmooth Weighted','Unweighted')


xlabel('Time (mins)')

set(gca,'FontSize',18)
set(gca,'FontName','Times New Roman')
hold off



