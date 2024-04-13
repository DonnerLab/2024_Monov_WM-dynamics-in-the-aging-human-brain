% Get data for demographic sample descriptions for Table 1 
% AGE | SEX | Diagnostic Subgroups MCI 
% Gina Monov, UKE, 2022 

clear all

% Load subject IDs
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/Final/behav_subj.mat']);
load(['/Users/ginamonov/Servers/mountpoint1/behavioral_analysis/Groups/yhc.mat']);
% Load data 
data_mci =  readtable('/Users/ginamonov/Servers/mountpoint1/psych_data/CRF_T1.csv'); 
data_yhc =  readtable('/Users/ginamonov/Servers/mountpoint1/sample_description/Datensatz_Entscheidungsstudie_0107.csv'); 
load('/Users/ginamonov/Servers/mountpoint1/psych_data/PCA_neuropsych_data/cerad_table_PCA.mat')

% Prepare MCI IDs in the table (without 'MCI_')
for t = 1:height(data_mci) 
    ID = []; 
    ID = data_mci.ID(t); 
    ID = ID {1,1}; 
    ID = ID(5:end); 
    ID = {ID}; 
    data_mci.ID(t) = ID; 
end 

hc = behav_hc_final; 
mci = behav_mci_final; 
unc = behav_cog_def_final; 

% Loop over all subgroups to get the corresponding data 

for s = 1:length(yhc)
   idx =  find(strcmp(data_yhc.ID,yhc{s}))
   age_yhc(s) = data_yhc.Alter(idx); 
   edu_yhc(s) = data_yhc.Bildung(idx);
   sex_yhc(s) = data_yhc.Geschlecht(idx);
end 


for s = 1:length(hc)
   idx =  find(strcmp(data_mci.ID,hc{s}))
  if strcmp(hc{s},'31') % manually add this for 31 since CRF data was not obtained (drop out) 
       age_hc(s)= 72; 
       sex_hc(s) = 2;     
  else

   age_hc(s) = data_mci.CRF_02(idx);

   edu_school_hc(s) = data_mci.CRF_14(idx);

   edu_prof_hc(s) = data_mci.CRF_14(idx);

   sex_hc(s) = data_mci.CRF_03(idx);
  end 
end 


for s = 1:length(mci)
   idx =  find(strcmp(data_mci.ID,mci{s})); 
   age_mci(s) = data_mci.CRF_02(idx);
   edu_school_mci(s) = data_mci.CRF_14(idx);
   edu_prof_mci(s) = data_mci.CRF_14(idx);
   sex_mci(s) = data_mci.CRF_03(idx);
   if strcmp(mci{s},'02') % manually add diagnosis for 2 because subject was excluded from PCA (not all tests were performed)
       diagnosis(s)= 2; 
   else 
   idx =  find(strcmp(cerad_data.ID,mci{s})); 
   diagnosis(s)= cerad_data.Diagnose(idx); 
   end 
end 

for s = 1:length(unc)
   idx =  find(strcmp(data_mci.ID,unc{s})); 
   age_unc(s) = data_mci.CRF_02(idx);
   edu_school_unc(s) = data_mci.CRF_14(idx);
   edu_prof_unc(s) = data_mci.CRF_14(idx);
   sex_unc(s) = data_mci.CRF_03(idx);
end 

% Switch numeric allocation of male and female in yhc to match the older subjects 
sex_yhc(sex_yhc==1) = 3; 
sex_yhc(sex_yhc==2) = 1;
sex_yhc(sex_yhc==3) = 2;

demo_table.YHC(1) = {[mean(age_yhc),std(age_yhc)]}; 
demo_table.YHC(2) = {length(sex_yhc(sex_yhc==2))./length(sex_yhc)*100}; 

demo_table.OHC(1) = {[mean(age_hc),std(age_hc)]}; 
demo_table.OHC(2) = {length(sex_hc(sex_hc==2))./(length(sex_hc))*100}; 

demo_table.MCI(1) = {[mean(age_mci),std(age_mci)]}; 
demo_table.MCI(2) = {length(sex_mci(sex_mci==2))./length(sex_mci)*100}; 

demo_table.UNC(1) = {[mean(age_unc),std(age_unc)]}; 
demo_table.UNC(2) = {length(sex_unc(sex_unc==2))./length(sex_unc)*100}; 
disp('multidomain-amnestic MCI, n =')
length(diagnosis(diagnosis==2))
disp('amnestic MCI, n =')
length(diagnosis(diagnosis==3))
