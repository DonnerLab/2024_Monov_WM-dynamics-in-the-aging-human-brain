% Compute group differences OHC|MCI for CERAD test battery subtests 
% --> Extended Data Table 1-1
% Permutation test: https://github.com/rudyvdbrink/Tools/tree/master/statistics
% Gina Monov, UKE, 2024

clear all

% Load CERAD data 
load('/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')

% Load subject IDs
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
behav_mci_final = behav_mci_final(~ismember(behav_mci_final,'02')); % Exclude subject that has not performed all subtests 
hc = behav_hc_final;
pat = behav_mci_final;
cog_def = behav_cog_def_final; 
subj = horzcat (hc,pat,cog_def);

% Load data for each subgroup 
ohc_counter = 0; 
mci_counter = 0; 
for l = 1:height(cerad_data) 
    
    if ismember(cerad_data.ID(l),hc)
       ohc_counter = ohc_counter +1; 
       WL_OHC(ohc_counter) = cerad_data.WL(l); 
       WLAB_OHC(ohc_counter) = cerad_data.WLAB(l); 
       VF_OHC(ohc_counter) = cerad_data.VF(l);
       BNT_OHC(ohc_counter) = cerad_data.BNT(l); 
       discrim_OHC(ohc_counter) = cerad_data.discrim(l); 
       VK_1_OHC(ohc_counter) = cerad_data.VK_1(l);
       VK_2_OHC(ohc_counter) = cerad_data.VK_2(l); 
       TMTA_OHC(ohc_counter) = cerad_data.TMTA(l);
       TMTB_OHC(ohc_counter) = cerad_data.TMTB(l); 
       PF_OHC(ohc_counter) = cerad_data.PF(l); 
    end 
    
    if ismember(cerad_data.ID(l),pat)
       mci_counter = mci_counter +1; 
       WL_MCI(mci_counter) = cerad_data.WL(l); 
       WLAB_MCI(mci_counter) = cerad_data.WLAB(l);
       VF_MCI(mci_counter) = cerad_data.VF(l); 
       BNT_MCI(mci_counter) = cerad_data.BNT(l);
       discrim_MCI(mci_counter) = cerad_data.discrim(l);
       VK_1_MCI(mci_counter) = cerad_data.VK_1(l);
       VK_2_MCI(mci_counter) = cerad_data.VK_2(l);
       TMTA_MCI(mci_counter) = cerad_data.TMTA(l);
       TMTB_MCI(mci_counter) = cerad_data.TMTB(l);
       PF_MCI(mci_counter) = cerad_data.PF(l);
    end 
    
end 

% Pull mean and standard deviation for each group 
WL.MCI_mean = mean(WL_MCI); 
WL.MCI_std = std(WL_MCI); 
WL.OHC_mean = mean(WL_OHC); 
WL.OHC_std = std(WL_OHC); 

WLAB.MCI_mean = mean(WLAB_MCI); 
WLAB.MCI_std = std(WLAB_MCI); 
WLAB.OHC_mean = mean(WLAB_OHC); 
WLAB.OHC_std = std(WLAB_OHC); 

BNT.MCI_mean = mean(BNT_MCI); 
BNT.MCI_std = std(BNT_MCI); 
BNT.OHC_mean = mean(BNT_OHC); 
BNT.OHC_std = std(BNT_OHC); 

VF.MCI_mean = mean(VF_MCI); 
VF.MCI_std = std(VF_MCI); 
VF.OHC_mean = mean(VF_OHC); 
VF.OHC_std = std(VF_OHC); 

discrim.MCI_mean = mean(discrim_MCI); 
discrim.MCI_std = std(discrim_MCI); 
discrim.OHC_mean = mean(discrim_OHC); 
discrim.OHC_std = std(discrim_OHC); 

VK_1.MCI_mean = mean(VK_1_MCI); 
VK_1.MCI_std = std(VK_1_MCI); 
VK_1.OHC_mean = mean(VK_1_OHC); 
VK_1.OHC_std = std(VK_1_OHC); 

VK_2.MCI_mean = mean(VK_2_MCI); 
VK_2.MCI_std = std(VK_2_MCI); 
VK_2.OHC_mean = mean(VK_2_OHC); 
VK_2.OHC_std = std(VK_2_OHC); 

TMTA.MCI_mean = mean(TMTA_MCI); 
TMTA.MCI_std = std(TMTA_MCI); 
TMTA.OHC_mean = mean(TMTA_OHC); 
TMTA.OHC_std = std(TMTA_OHC); 

TMTB.MCI_mean = mean(TMTB_MCI); 
TMTB.MCI_std = std(TMTB_MCI); 
TMTB.OHC_mean = mean(TMTB_OHC); 
TMTB.OHC_std = std(TMTB_OHC); 

PF.MCI_mean = mean(PF_MCI); 
PF.MCI_std = std(PF_MCI); 
PF.OHC_mean = mean(PF_OHC); 
PF.OHC_std = std(PF_OHC); 

% Compare subgroup performances 
 [h, p.WL, diff, diff_null]=permtest2(WL_OHC,WL_MCI,10000);
 [h, p.WL_recall, diff, diff_null]=permtest2(WLAB_OHC,WLAB_MCI,10000);
 [h, p.BNT, diff, diff_null]=permtest2(BNT_OHC,BNT_MCI,10000);
 [h, p.VF, diff, diff_null]=permtest2(VF_OHC,VF_MCI,10000);
 [h, p.discrim, diff, diff_null]=permtest2(discrim_OHC,discrim_MCI,10000);
 [h, p.CP, diff, diff_null]=permtest2(VK_1_OHC,VK_1_MCI,10000);
 [h, p.CP_recall, diff, diff_null]=permtest2(VK_2_OHC,VK_2_MCI,10000);
 [h, p.TMTA, diff, diff_null]=permtest2(TMTA_OHC,TMTA_MCI,10000);
 [h, p.TMTB, diff, diff_null]=permtest2(TMTB_OHC,TMTB_MCI,10000);
 [h, p.PF, diff, diff_null]=permtest2(PF_OHC,PF_MCI,10000);
 
 