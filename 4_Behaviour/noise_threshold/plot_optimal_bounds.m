% Plot accuracy optimizing threshold as a function of memory noise in synthetic data set 
% --> Extended Data Figure 3-2A
% Gina Monov, UKE, 2023

clear all
close all

load(['/Users/ginamonov/Servers/mountpoint1/Modelling/Dec_Rule_Modelling/optimal_criterion/','optimal_bounds.mat'])
load('/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat')

for l = 1:length(optimal_bounds) 
    noise(l) = optimal_bounds(l).noise;
    criterion(l)=optimal_bounds(l).optimal; 
    
end 

h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 5; % figure width
fig_h = 5; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

subplot(1,1,1), hold on 
axis square 
ylabel('Optimal threshold \delta')
xlabel('Noise \sigma_{mem}')

[rho,pval]=corr([noise]',[criterion]','type','Pearson');       
p = polyfit([noise]',[criterion]',1); 
f = polyval(p,[noise]'); 
scatter(noise,criterion,20,'MarkerFaceColor',colors.white,'MarkerEdgeColor',colors.grey,'MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5), hold on 
plot(noise,f,'Color',colors.black,'LineWidth',2), hold on
 
if pval < 10^-4
  txt = {['r=' num2str(round(rho,4))]; ['p<10^{-4}']};
  text(2,11,txt,'Color',colors.black,'fontsize',7,'fontweight','bold'), hold on 
else
  text(2,11,txt,'Color',colors.grey,'fontsize',7), hold on 
end 
set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',7,'TickDir','out','box','off')

cd '/Users/ginamonov/Servers/mountpoint1/final_figures/'   
savefig(figure(1),['optimal_criterion.fig']) % Extended Data Figure 3-2A