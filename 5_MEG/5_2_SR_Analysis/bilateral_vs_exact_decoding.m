% Plot bilateral and exact stimulus onset decoding for 1 and 3s in dorsal
% visual cortex 
% cluster-based permutation test: Edden M. Gerber (2023). permutest, MATLAB Central File Exchange. https://www.mathworks.com/matlabcentral/fileexchange/71737-permutest
% --> Extended Data Fig. 6-1A
% Gina Monov, UKE 2024 

clear all 
close all
addpath '/Users/ginamonov/Servers/mountpoint1/functions/'
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/meg_subj.mat']);
load('/Users/ginamonov/Servers/mountpoint1/functions/colors/colors.mat')

mci = meg_mci_final; 
hc = meg_hc_final; 
cog_def = meg_cog_def_final;
subj = horzcat(mci,hc,cog_def);

delays = [1,3]; 

% Initialize figure 
h =  findobj('type','figure');
cfig = length(h)+1;
fig_w = 14; % figure width
fig_h = 6; % figure height
f1=figure; set(cfig,'units','centimeters','pos', [0 0 fig_w fig_h],'Color',[1 1 1])
set(gcf,'renderer','Painters'),

for d = 1:length(delays) 

for s = 1:length(subj) 
 
 
    exact_loc = readtable(['/Users/ginamonov/Servers/mountpoint1/meg_analysis/Decoding/decode_5_35_coarse/' subj{s}, '_Dorsal_visual_', num2str(delays(d)),'_5_35.csv']); 
    bi_loc = readtable(['/Users/ginamonov/Servers/mountpoint1/meg_analysis/Decoding/bi_decode_5_35_coarse/bi_' subj{s}, '_Dorsal_visual_', num2str(delays(d)),'_5_35.csv']);
   
    times = bi_loc.latency; 
   
             v = 1:length(times); 
             oddidx = @(v) v(1:2:end);           % Addressing Odd-Indexed Element since each time point was accidentally saved twice
             y1 = oddidx(v); 
             exact(s,:) = exact_loc.test_correlation(y1); 
             bi(s,:) = bi_loc.test_score(y1); 
             times = times(y1); 

    clear exact_loc bi_loc 
     
end 
 

% Perform cluster-based permutation test against zero for decoding time course 

data4stats = []; 
data4stats = exact'; %Flip around to bring it into the right format 
null4stats = zeros(size(data4stats)); 
[clusters,p_values,t_sums,permutation_distribution] = permutest(data4stats,null4stats,true,0.05,10000,false);

yyaxis left %plot correlation coefficients on the left
   
    % Plotting 
    
subplot(1,2,d), hold on 
    s1=plot(times,mean(exact,1),'Color',colors.sky), hold on 
    shadedErrorBar(times,mean(exact,1),std(exact,[],1)./sqrt(size(exact,1)),{'Color',colors.sky},1)
         
         % add statistics to plots 
         for ccc = 1:length(clusters) 
             if p_values(ccc) < 0.05
             plot(times(clusters{1,ccc}),zeros(size(clusters{1,ccc}))-0.01,'Color',colors.sky,'LineWidth',1.5), hold on 
             end
         end 
         
 if d == 1        
    set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.5 1 1.5],'YTick',[0 0.1 0.2 0.3 0.4])
    set(gca,'XTickLabel',{'0' '0.5' '1' '1.5'},'YTickLabel',{ '0'  '0.1' '0.2'  '0.3' '0.4'})
 elseif d == 2
    set(gca,'FontName','Helvetica','LineWidth',0.5,'FontSize',8,'TickDir','out','box','off','XTick',[0 0.5 1 2 3 3.5],'YTick',[0 0.1 0.2 0.3 0.4])
    set(gca,'XTickLabel',{'0' '0.5' '1' '2' '3' '3.5'},'YTickLabel',{ '0'  '0.1' '0.2'  '0.3' '0.4'})    
 end 
 
onsets = [0 0.5]; %mem on and mem off
xline(onsets(1),'--','Color',colors.black,'LineWidth',1), hold on 
xline(onsets(2),'--','Color',colors.black,'LineWidth',1), hold on 
yline(0,'-','Color',colors.black), hold on % Reference line around zero
ylim([-0.05  0.4])
ylabel('Correlation coefficient','Fontsize',8)         
                   
data4stats = []; 
data4stats = bi'; %Flip around to bring it into the right format 
null4stats = zeros(size(data4stats))+0.5; 
[clusters,p_values,t_sums,permutation_distribution] = permutest(data4stats,null4stats,true,0.05,10000,false);

         % add statistics to plots 
         for ccc = 1:length(clusters) 
             if p_values(ccc) < 0.05
             plot(times(clusters{1,ccc}),zeros(size(clusters{1,ccc}))-0.025,'Color',colors.tiger,'LineWidth',1.5), hold on 
             end
         end 

 % Plotting 
  yyaxis right % plot auc values to the right 

  xlim([-0.1 delays(d)+0.5]); ylim([0.5-(0.05.*(0.25./0.4))  0.75]) 
  s2=plot(times,mean(bi,1),'Color',colors.tiger), hold on 
  shadedErrorBar(times,mean(bi,1),std(bi,[],1)./sqrt(size(bi,1)),{'Color',colors.tiger},1), hold on 
         
xlabel('Time from stimulus onset (s)','Fontsize',8)
ylabel('Area under the ROC-curve','Fontsize',8,'Color',colors.black)    
ax = gca; 
ax.YColor = 'k';
if d == 1         
sgtitle('Dorsal visual','Fontsize',10,'Color', colors.black)
legend([s1,s2],{'exact sample location','sample hemifield'},'fontsize',8)
end 
clear bi exact
hold on 
end 