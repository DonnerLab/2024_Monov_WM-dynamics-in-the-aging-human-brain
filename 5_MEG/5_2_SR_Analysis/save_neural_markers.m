% Compute and save mean decoding precision during delay as neural marker for each subject and
% Glasser parcel (N = 180)
% Gina Monov, UKE, 2023

clear all
close all
load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
allsubj={}; 
for f = 1:length(datasets_overview) 
    if datasets_overview(f).sess == 1
    allsubj = horzcat(allsubj,datasets_overview(f).subj);
    end 
end 
load('/mnt/homes/home028/gmonov/meg_analysis/all_glasser_parcels.mat')
new_ROIs={};
for g = 1:length(ROIs) 
    new_ROIs{g} = replace(ROIs{g},'_','-'); 
end 
ROIs = new_ROIs; 
delay = '1'; 
 

for s = 1:length(allsubj)

   for l = 1:length(ROIs) 
     try  
     dec_results = readtable(['/mnt/homes/home028/gmonov/meg_analysis/Decoding/decode_5_35_complete_glasser/' allsubj{s} '_' ROIs{l} '_' delay '_5_35.csv']);
     if ~isnan(dec_results.test_correlation(1))
         times = dec_results.latency; 
         corr = dec_results.test_correlation; 
         dec_prec = mean(corr(times>= 0.5 & times <1.5)); 
         save(['/mnt/homes/home028/gmonov/meg_analysis/neural_markers/all_parcels/dec_prec/',allsubj{s}, '_', ROIs{l}, '.mat'], 'dec_prec');
     else disp([allsubj{s}, ', ', ROIs{l}, ' nan decoding results'])
     end 
     catch disp([allsubj{s}, ', ', ROIs{l}, ' no decoding results'])
     end 
   end 
end 