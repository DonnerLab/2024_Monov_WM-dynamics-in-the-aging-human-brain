% Generates behaviour from Wiener diffusion model of
% delayed-match-to-sample working memory task with sigmoidal decision rule. 
% Model instances are determined to approximate their value within the
% fitted behaviour. 
% pm: model parameters ([sigmaM sigmaD bound theta lambda])
% delays: vector of delay period durations (real task: [1 3 9])
% tpd: number of trials per delay (real task: 21 per block)
% p_match: proportion of trials on which sample matches probe (real task: 0.33)
% p_near: proportion of trials on which probe is presented at nearest location (real task: 0.33)
% All inputs currently hardcoded within the function, need to be commented
% for usage within other scripts.
% Peter Murphy, Trinity College Dublin & Gina Monov, UKE, 2022

function [dat] = gen_data_wd_DecRule(pm,delays,tpd,nblocks,p_match,p_near,n_datasets) 


modeltype = {'sM+bound','sM+bound+theta','sM+bound+lambda','sM+sD+bound','sM+bound+theta+lambda','sM+sD+bound+theta+lambda'}; %add all parameter settings
load(['/mnt/homes/home028/gmonov/behavioral_analysis/Groups/yhc.mat']);
loadpath_mci = '/mnt/homes/home028/gmonov/Modelling/Dec_Rule/'; 
loadpath_scz = '/mnt/homes/home028/gmonov/SCZ/WM_modeling/Dec_Rule/'; 

%Pull subject IDs
load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
allsubj={};
for f = 1:length(datasets_overview) 
    if datasets_overview(f).sess == 1
    allsubj = horzcat(allsubj,datasets_overview(f).subj);
    end 
end 
subj = horzcat(allsubj,yhc); %all subjects in the mci study including the excluded ones
try
load(['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/gen_datasets/model_instances/median_fits.mat'])
catch
gen_pm = [nan,nan,nan,nan,nan]; % Pull model instances for parameter recovery


for t = 1:length(modeltype) 
    
    sigmaM_y = []; %Initialize all possible free parameters
    sigmaD_y = [];
    bound_y = [];
    theta_y = [];
    lambda_y = [];
    
    for s = 1:length(subj)
        
    if sum(ismember(subj{s},allsubj))==1
    fullpath = ['/mnt/homes/home028/gmonov/Modelling/Dec_Rule/',modeltype{t}];
    elseif sum(ismember(subj{s},yhc))==1
    fullpath = ['/mnt/homes/home028/gmonov/SCZ/WM_modeling/Dec_Rule/',modeltype{t}];
    end
        load([fullpath,filesep,subj{s},'.mat']);
 
      sigmaM = []; %Initialize all possible free parameters
      sigmaD = [];
      bound = [];
      theta = [];
      lambda = [];
    
       for r=1:length(pm)
  
         if isempty(sigmaM) & ismember('sM',modeltype{t})
         sigmaM_y(s) = pm(r);
         sigmaM = pm(r); 
         elseif isempty(sigmaD) & ismember('sD',modeltype{t})
         sigmaD_y(s) = pm(r);
         sigmaD = pm(r); 
         elseif isempty(bound) & ismember('bound',modeltype{t})
         bound_y(s) = pm(r);
         bound = pm(r); 
         elseif isempty(theta) & ismember('theta',modeltype{t})
         theta_y(s) = pm(r);
         theta = pm(r); 
         elseif isempty(lambda) & ismember('lambda',modeltype{t})
         lambda_y(s) = pm(r);
         lambda = pm(r); 
         end 
     
      end 
    
    end 
    
    %Get median of all subjects to work with for parameter recovery 
     if isnan(gen_pm(1))
       gen_pm(1) = median(sigmaM_y); 
     end
     if ~isempty(sigmaD_y) && isnan(gen_pm(2))
       gen_pm(2) = median(sigmaD_y);
     end
     if isnan(gen_pm(3))
        gen_pm(3) = median(bound_y);
     end
     if ~isempty(theta_y) && isnan(gen_pm(4))
      gen_pm(4) = median(theta_y);
     end 
     if ~isempty(lambda_y) && isnan(gen_pm(5))
     gen_pm(5) = median(lambda_y);
     end 
end 

save(['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/gen_datasets/model_instances/median_fits.mat'] , 'gen_pm');

end

clearvars -except gen_pm modeltype 


for model = 1:length(modeltype)
    
   modelspec = modeltype{model}; 

   model_log = zeros(1,5); % logical of free parameters in specified model


% Get logical which free parameters exist in this model
if ismember('sM',modelspec)
model_log(1) = 1; 
end 
if ismember('sD',modelspec)
model_log(2) = 1; 
end
if ismember('bound',modelspec)
model_log(3) = 1; 
end
if ismember('theta',modelspec)
model_log(4) = 1; 
end
if ismember('lambda',modelspec)
model_log(5) = 1; 
end
   pm = gen_pm; % renaming 
% Get parameters
   if model_log(1) == 1
    sigmaM = pm(1);
   else sigmaM = [];
   end 
   if model_log(2) == 1
    sigmaD = pm(2);
   else sigmaD = [];
   end 
   if model_log(3) == 1
    bound = pm(3);
   else bound = [];
   end 
   if model_log(4) == 1
    theta = pm(4);
   else theta = [];
   end 
   if model_log(5) == 1
    lambda = pm(5);
   else lambda = [];
   end 
    %only feed in the ones specified in modelspec
 
    n_datasets = 100; % How many datasets to generate?

for sets = 1:n_datasets
    


delays = [1 3 9];
tpd = 21; 
p_match = 0.33; 
p_near = 0.33; 
nblocks = 3; 

tpd = tpd*nblocks; %increase trials per delay corresponding to desired number of blocks 




% Exact locations (polar angles) used in task (only used here except for calculating sample-probe distances)
nlocs = 12;
angles = [0:(180/(nlocs+1)):180];  % possible stimulus polar angles (adding 2 so that most extreme memoranda can still be flanked either side by a near non-match)
anglesU = angles(2:end-1);

% Seed random number generator
rng('default')
rng('shuffle')

%%% Generate trial info using same process as we use in task script
% MATCH trials
matchS = nan(round(tpd*p_match),length(delays)); matchD = [];
for d = 1:length(delays)
    if round(tpd*p_match) == length(anglesU)
        matchS(:,d) = anglesU';
    elseif round(tpd*p_match) > length(anglesU)
        tempmatchS = repmat(anglesU',floor(round(tpd*p_match)/length(anglesU)),1);
        tempmatchS = [tempmatchS; randsample(anglesU,round(mod(round(tpd*p_match),length(anglesU))))'];
        matchS(:,d) = tempmatchS;
    else matchS(:,d) = randsample(anglesU,round(tpd*p_match))';
    end
    matchD = [matchD; ones(size(matchS,1),1).*delays(d)];
end
matchS = reshape(matchS,length(matchD),1);

% Construct matrix of NEAR-NO-MATCH sample stimuli (ensuring at least one of each sample stimulus per delay if possible)
nearS = nan(round(tpd*p_near),length(delays)); nearD = [];
for d = 1:length(delays)
    if round(tpd*p_near) == length(anglesU)
        nearS(:,d) = anglesU';
    elseif round(tpd*p_near) > length(anglesU)
        tempnearS = repmat(anglesU',floor(round(tpd*p_near)/length(anglesU)),1);
        tempnearS = [tempnearS; randsample(anglesU,round(mod(round(tpd*p_near),length(anglesU))))'];
        nearS(:,d) = tempnearS;
    else nearS(:,d) = randsample(anglesU,round(tpd*p_near))';
    end
    nearD = [nearD; ones(size(nearS,1),1).*delays(d)];
end
nearS = reshape(nearS,length(nearD),1);

% Construct matrix of FAR-NO-MATCH sample stimuli (ensuring at least one of each sample stimulus per delay if possible)
p_far = 1-p_match-p_near;
farS = nan(round(tpd*p_far),length(delays)); farD = [];
for d = 1:length(delays)
    if round(tpd*p_far) == length(anglesU)
        farS(:,d) = anglesU';
    elseif round(tpd*p_far) > length(anglesU)
        tempfarS = repmat(anglesU',floor(round(tpd*p_far)/length(anglesU)),1);
        tempfarS = [tempfarS; randsample(anglesU,round(mod(round(tpd*p_far),length(anglesU))))'];
        farS(:,d) = tempfarS;
    else farS(:,d) = randsample(anglesU,round(tpd*p_far))';
    end
    farD = [farD; ones(size(farS,1),1).*delays(d)];
end
farS = reshape(farS,length(farD),1);

% Construct matrices of PROBE stimuli
matchP = matchS;

nearP = nearS;
for i = 1:size(nearP,1)
    c_angles = [find(angles==nearS(i))-1 find(angles==nearS(i))+1];
    nearP(i) = angles(randsample(c_angles,1));
end

farP = farS;
for i = 1:size(farP,1)
    no_angle = find(angles==farS(i))-1:find(angles==farS(i))+1;
    c_angles = find(~ismember(1:length(angles),no_angle));
    farP(i) = angles(randsample(c_angles,1));
end

% Construct full matrix of stimulus types [SAMPLE ANGLE; PROBE ANGLE; DELAY; TRIAL-TYPE] (trial-types: 1= match, 2=near no-match, 3=far no-match)
stimIn = [[matchS matchP matchD ones(size(matchP,1),1)]; [nearS nearP nearD ones(size(nearP,1),1).*2]; [farS farP farD ones(size(farP,1),1).*3]];

% Randomize order
shuforder = randperm(size(stimIn,1));
stimIn = stimIn(shuforder,:);

% Calculate vectors of single-trial delay durations and M-P distances
delays = stimIn(:,3);
deltas = round(abs(stimIn(:,2)-stimIn(:,1)),1); % rounding because otherwise there are small rounding errors which make indexing annoying


%%% Generate single-trial choice probabilities and choices
for t = 1:size(stimIn,1)
    cp(t,1) = run_wd_DecRule(sigmaM,sigmaD,bound,theta,delays(t),deltas(t));  % CP given memory noise and decision rule
    if ~isempty(lambda)
    cp(t,1) = cp(t).*(1-(1-exp(-lambda*delays(t)))) + 0.5.*(1-exp(-lambda*delays(t)));  % add memory lapse component
    end
end
resp = binornd(ones(size(cp)),cp);   % draw binary choices from choice probabilities


%%% Collate output into format required for model fitting
dat = [delays deltas resp];


save(['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/gen_datasets/',modelspec, '/', num2str(sets) '.mat'] , 'dat');

clearvars -except modeltype model sets gen_pm sigmaM sigmaD bound theta lambda modelspec
end 
end