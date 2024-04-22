# 2024_Monov_WM-dynamics-in-the-aging-human-brain

📄 This project contains codes from the publication: Monov, G., Stein, H., Klock, L., Gallinat, J., Kühn, S., Lincoln, T., Krkovic, K., Murphy, P.†, Donner, T.†, 2024. Linking cognitive integrity to working memory dynamics in the aging human brain (in press). Journal of Neuroscience. [† contributed equally]. [Preprint available here](https://www.biorxiv.org/content/10.1101/2023.08.18.553840v1).

It contains several folders with code for preprocessing and analyzing behavioral, MEG and eye-tracking data. Below is a short description of the contents within each subfolder and pointers to specific analyses and figures from the publication. Detailed descriptions can be found in the beginning of each of the scripts. Please [get in touch](mailto:monov.gina@gmail.com) if you have any questions and/or requests!
_________________________________________________________________

## 1_Neuropsych_and_sample_description/ 👵👴
Analyses of participant demographic data, cognitive performance scores and data visualization. 

•	`demographic_sample_description.m`: *Table 1.*

•	`group_comparisons_cerad_subtests.m`: *Extended Data Table 1-1.*

•	`run_final_pca.m`: *Figure 1B, Extended Data Figure 1-1A-B.*
_________________________________________________________________

## 2_WM_task/ 🔘👆
Working memory task and training protocol as performed in the MEG (‘Dataset 1’ in the publication), both coded in [Psychtoolbox-3 for Matlab](http://psychtoolbox.org/).

•	`WM_dm2s_block.m`: *Main script task block.*

•	`WM_dm2s_training.m`: *Main script training.*
________________________________________________________________
## 3_WM_model/ 👩🏻‍💻👨🏻‍💻
Fits algorithmic model of working memory to behavioral data from the task and generates synthetic data for model validation and comparison. 

### 3_1_Model_fitting/
•	`FIT_wm_diffusion_DecRule.m`: *Highest-level model fitting script.*
### 3_2_validation_comparison/
•	`Dec_Rule_model_comparsion.m`: *Extended Data Fig. 3-1C.*

•	`gen_data_wd_DecRule.m`: *Creates synthetic data for parameter recovery.*

•	`parameter_recovery_plots.m`:  *Extended Data Fig. 3-1A-B.*

•	`R2_behavioral_model.m`: *Extended Data Fig. 3-1E.*

•	`save_miss_fa_H.m` & `SDT_model_params_corr.m`: *Extended Data Fig. 3-1F.*

•	`visualizations_Dec_rule.m`: *Figure 3A, Extended Data Fig. 3-1D.*
_________________________________________________________________

## 4_Behaviour/ 📊
Analyses of behavioral data (Model-free and model-based). 

•	`behav_plots4paper.m`: *Figure 2, Figure 3B, Extended Data Fig. 2; computes and saves lapse+noise parameters, Mixed ANOVA analyses.*

•	`bf_anova_model_params.m`: *Further analyses on group differences of fitted model parameters.*

•	`check_behavior.m`: *Extended Data Fig. 1-1C; group-wise comparison of task accuracy.*

### noise_threshold/
•	`NEW_optimal_criterions.m` & `find_optimal_bound.m`: *Fit optimal threshold parameters for different levels of memory noise.*

•	`plot_optimal_bounds.m`: *Extended Data Figure 3-2A.*

•	`thres_noise_acc_corr.m`: *Extended Data Figure 3-2B.*

### PC1_corr/
•	`Figure_4_4paper.m`: *Figure 4, Extended Data Fig. 4-1 A-B.*

•	`check_corrs_with_cerad_total.m`: *Extended Data Fig. 4-1C.*

_________________________________________________________________

## 5_MEG/ 🧠
Preprocessing and source reconstruction/analysis of MEG data (Dataset 1). 

### 5_1_Preprocessing/
•	`dataset_info/…`: *Preparatory scripts for creating master file containing information for each dataset that can be pulled later.*

•	`meg_preproc/…`: *Further preprocessing of MEG data proceeded according to a [standardized pipeline developed in our laboratory](https://github.com/DonnerLab/meg-preproc). Current pipeline was lightly altered to conform to the goals of the present study.*

•	`prep_SR/…`: *Final steps to prepare data for source reconstruction performed in MNE Python.*

•	`sanity_checks.m` & `datasets_checks.m`: *Final checks ensuring that behavioral and MEG data are aligned, adding further information to master file and saving final master file.*

### 5_2_SR_Analysis/
Scripts for source reconstruction and analyses of the MEG data. 

•	`pymeg/…`: *Code for source reconstruction and decoding analysis as previously developed in [Wilming et al. 2020](https://www.nature.com/articles/s41467-020-18826-6) & [Murphy et al. 2021](https://www.nature.com/articles/s41593-021-00839-z). Uploaded codes have been specifically edited to conform the purpose of our study.*

•	`analyze_temp_gen.m`: *Figure 8C-E, Extended Data Figure 8-1.*

•	`bilateral_vs_exact_decoding.m`: *Extended Data Figure 6-1A.*

•	`brainplots4paper.m`: *Extended Data Figure 5-1B, 6-1B, 7-1A-B.* 

•	`brainplots4paper22.m`: *Figure 5B, Figure 6B.*

•	`Decode_plots4paper.m`: *Figure 6A.* 

•	`plot_temp_gen.m`: *Figure 8A-B.*

•	`roi_depiction.m` & `TFR_plots4paper.m`: *Figure 5A.* 

•	`TFRs_FigureS5_4_paper.m`: *Extended Data Figure 5-1A.* 

•	`WHOLE_CORTEX_multiple_linear_regr_neural_markers.m`: *Figure 7.*

_________________________________________________________________

## 6_Eye_tracking/ 👀
Preprocessing and analysis of eye tracking data (Datasets 1 & 2). 

### 6_1_Preprocessing/
Contains separate subfolders for Datasets 1 and 2. The main scripts are denoted as `S1_*.m`, `S2_*.m`, `S3_*.m`, where the numbers indicate the order in which they are executed. 

•	`eye_prep_trial_analysis.m`: *Prepares data of both datasets for decoding analysis and checks alignment with behavioral data.*

### 6_2_Analysis/
•	`eye_decoding/…`: *Decoding of sample stimulus location on gaze position.*

•	`analyze_eye_dec_1000Hz.m`: *Extended Data Figure 9-1.* 

•	`analyze_eye_dec.m`: *Figure 9.*

•	`corr_eye_meg_dec.m`: *Correlates MEG decoding precision for each Glasser group 
with gaze decoding precision.* 

_________________________________________________________________
Please note that some scripts make use of functions not created by our team and that are available for download elsewhere. These are referred to at the beginning of our scripts where appropriate. Below is a table listing these functions, including a description of their use for the publication and an external link to each of them. 

| Reference  | Used in |
| ------------- | ------------- |
| [Bart Krekelberg (2023). bayesFactor. GitHub.](https://github.com/klabhub/bayesFactor)  | Computation of bayes factor of group differences for fitted threshold parameters  |
| [Bechtold, Bastian (2016). Violin Plots for Matlab. DOI: 10.5281/zenodo.4559847](https://github.com/bastibe/Violinplot-Matlab)  | Violin plots of working memory accuracy in Extended Data Figure 1-1C  |
| [Birge B (2003). PSOt - a particle swarm optimization toolbox for use with Matlab. In: Proceedings of the 2003 IEEE Swarm Intelligence Symposium. SIS’03 (Cat. No.03EX706), pp 182–186. Indianapolis, IN, USA: IEEE.](http://ieeexplore.ieee.org/document/1202265/)  | Particle Swarm Optimization Algorithm for behavioral model fitting |
| [D’Errico J (2012). fminsearchbnd, fminsearchcon. MATLAB Central File Exchange.](https://www.mathworks.com/matlabcentral/fileexchange/8277-fminsearchbnd-fminsearchcon)  | Simplex minimization for fitting optimal threshold for noise levels  |
| [Droste Effect – Brewermap. Github.](https://github.com/DrosteEffect/BrewerMap)  | Colormaps for brain maps and temporal generalization matrices  |
| [Edden M. Gerber (2023). permutest. MATLAB Central File Exchange.](https://www.mathworks.com/matlabcentral/fileexchange/71737-permutest)   | Cluster-based permutation test of decoding precision time courses and temporal generalization matrices |
| [Laurent Caplette (2019). Simple RM/Mixed ANOVA for any design. MATLAB Central File Exchange.](https://www.mathworks.com/matlabcentral/fileexchange/64980-simple-rm-mixed-anova-for-any-design)  | ANOVA analyses of model-free behavioral analyses |
| [Philip Yip (2020). colorpalette. MATLAB Central File Exchange.](https://www.mathworks.com/matlabcentral/fileexchange/68969-colorpalette)   | Color definition in most of the analysis codes |
| [RL van den Brink. Tools/statistics. Github.](https://github.com/rudyvdbrink/Tools/tree/master/statistics)  | Non-parametric permutation tests for within-subjects comparisons against zero and between-subjects group-wise comparisons|
| [Rob Campbell (2019). raacampbell/shadedErrorBar. Github.](https://github.com/raacampbell/shadedErrorBar)  | Shaded error bars |
| [Scott Lowe (2022). superbar. MATLAB Central File Exchange.](https://www.mathworks.com/matlabcentral/fileexchange/57499-superbar)  | Bar plots and corresponding error bars  |
