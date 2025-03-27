%%%%
%%%% Analysis of fitting results for the STAR synthetic dataset.

clear all;
clc;
close all
global par p


load STAR_Nominal_Parameter.mat;

True_p = p;

%%% Load Results %%%
load Simulated_Data_Jan22_2025.mat;
True = Simu_Traje';

clear Simu_Traje p

Simu_Traje = True;
clear True


% Initial Values
par.PYtot = 0.5*10^(-9);  %%% 0.5nM STAR targeting GFP plasmid

par.Ps = 8*10^(-9); %%% STAR 8nM


load STAR_Non_Smooth_All_Results.mat;
load STAR_Weight_Smooth_All_Results.mat;

%%%% Weighted -- GREEN; Weighted with smoothing -- RED; Unweighted -- blue


colors = [0 0.4470 0.7410; 0.4660 0.6740 0.1880; 0.8500 0.3250 0.0980];

Hcolors = [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880; 0 0.4470 0.7410];

% % yellow -- [0.9290 0.6940 0.1250];
% % red -- [0.8500 0.3250 0.0980]
% % green -- [0.4660 0.6740 0.1880]
%%% blue -- [0 0.4470 0.7410]


figure
for ij = 1:8
    subplot(2,4,ij)
    boxplot([Unweighted_p(:,ij) Weighted_p(:,ij) Weighted_Smooth_p(:,ij)],[1 2 3],...
        'BoxStyle','outline','Colors',colors,'Symbol','o','OutlierSize',6);
    Hbox = gca();
    set(Hbox.Children.Children,'LineWidth',2)
    set(Hbox,'XTickLabel',{'UnW';'Wt';'WtS'})
   
    h = findobj(gca,'Tag','Box');
    for j=1:length(h)
        patch(get(h(j),'XData'),get(h(j),'YData'),Hcolors(j,:),'FaceAlpha',.5);
    end

    hold on
    yline(True_p(ij),'-.','LineWidth',2);

    if ij == 1
        ylabel('\alpha_s')
    elseif ij == 2
        ylabel('\delta_s')
    elseif ij == 3
        ylabel('\delta_m')
    elseif ij == 4
        ylabel('\beta_s')
    elseif ij == 5
        ylabel('\alpha_m')
    elseif ij == 6
        ylabel('KI')
    elseif ij == 7
        ylabel('KE')
    else
        ylabel('\alpha_{gm}')
    end

    set(gca,'FontSize',20)
    set(gca,'FontName','Times New Roman')
    hold off

end


WeightMatrix = zeros(3,49);
UnweightMatrix = zeros(3,49);
WeightSmoothMatrix = zeros(3,49);
TrueMatrix = zeros(3,49);

for i = 1:49
    WeightMatrix(1,i) = min(Weighted_Simu(i,:));  %%% 1 is for 4nm Weighted
    WeightMatrix(2,i) = mean(Weighted_Simu(i,:));
    WeightMatrix(3,i) = max(Weighted_Simu(i,:));
    
    UnweightMatrix(1,i) = min(Unweighted_Simu(i,:));  %%% 1 is for 4nm Weighted
    UnweightMatrix(2,i) = mean(Unweighted_Simu(i,:));
    UnweightMatrix(3,i) = max(Unweighted_Simu(i,:));
    
    WeightSmoothMatrix(1,i) = min(Weighted_Smooth_Simu(i,:));  %%% 1 is for 4nm Weighted
    WeightSmoothMatrix(2,i) = mean(Weighted_Smooth_Simu(i,:));
    WeightSmoothMatrix(3,i) = max(Weighted_Smooth_Simu(i,:));
    
    
    TrueMatrix(1,i) = min(Simu_Traje(i,:));  %%% 1 is for 4nm Weighted
    TrueMatrix(2,i) = mean(Simu_Traje(i,:));
    TrueMatrix(3,i) = max(Simu_Traje(i,:));
end



figure

subplot(1,3,1)
plot(Simu_t,UnweightMatrix(2,:),'Color',colors(1,:),'LineWidth',1.5)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [UnweightMatrix(1,:),fliplr(UnweightMatrix(3,:))];
fill(tx,py,colors(1,:),'FaceAlpha',0.3,'EdgeColor','none');
hold on

plot(Simu_t,TrueMatrix(2,:),'k','LineWidth',1.5)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [TrueMatrix(1,:),fliplr(TrueMatrix(3,:))];
fill(tx,py,'k','FaceAlpha',0.1,'EdgeColor','none');
hold on

xlim([0 240])
ylim([0 0.105])
ylabel('EGFP Conc. (\muM)')
%xticks([])
xlabel('Time (mins)')
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
hold off


subplot(1,3,2)
plot(Simu_t,WeightMatrix(2,:),'Color',colors(2,:),'LineWidth',1.5)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [WeightMatrix(1,:),fliplr(WeightMatrix(3,:))];
fill(tx,py,colors(2,:),'FaceAlpha',0.3,'EdgeColor','none');
hold on

plot(Simu_t,TrueMatrix(2,:),'k','LineWidth',1.5)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [TrueMatrix(1,:),fliplr(TrueMatrix(3,:))];
fill(tx,py,'k','FaceAlpha',0.1,'EdgeColor','none');
hold on

xlim([0 240])
ylim([0 0.105])
%ylabel('EGFP Conc. (\muM)')
xlabel('Time (mins)')
%xticks([])
set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
hold off




subplot(1,3,3)
plot(Simu_t,WeightSmoothMatrix(2,:),'Color',colors(3,:),'LineWidth',1.5)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [WeightSmoothMatrix(1,:),fliplr(WeightSmoothMatrix(3,:))];
fill(tx,py,colors(3,:),'FaceAlpha',0.3,'EdgeColor','none');
hold on

plot(Simu_t,TrueMatrix(2,:),'k','LineWidth',1.5)
hold on
tx = [Simu_t',fliplr(Simu_t')];
py = [TrueMatrix(1,:),fliplr(TrueMatrix(3,:))];
fill(tx,py,'k','FaceAlpha',0.1,'EdgeColor','none');
hold on

xlim([0 240])
ylim([0 0.105])
% ylabel('EGFP Conc. (\muM) - 8nM')
xlabel('Time (mins)')

set(gca,'FontSize',20)
set(gca,'FontName','Times New Roman')
hold off




