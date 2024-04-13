% Generates behaviour from Wiener diffusion model of
% delayed-match-to-sample working memory task with sigmoidal decision rule. 
% Threshold that produces highest accuracy for a given level of memory noise is found through simplex
% minimization 
% pm: optimal threshold for a given noise level 
% Model fitting algorithm: D'Errico J (2012) fminsearchbnd, fminsearchcon. MATLAB Central File Exchange Retrieved February 6, 2012.
% Gina Monov, UKE, 2023


addpath('/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/')
addpath('/mnt/homes/home028/gmonov/functions/')
addpath('/mnt/homes/home028/gmonov/functions/FMINSEARCHBND')

modeltype = {'sM+bound'}; % for this purpose only version with sM and bound

lambda = [];
theta = [];
sigmaD = [];

noise_levels = linspace(0.4,11,100); % Choose noise levels to test (app. range that was observed in subjects' model fits) 

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

% Generate trial info using same process as we use in task script
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


for sets = 1:length(noise_levels)
 
    sigmaM = noise_levels(sets);
    
          % Fitting initialization [r0 r1 beta kappa theta]
          inits = [10];  % starting point
          lb = [0];    % lower bound
          ub = [50];   % upper bound
          
% Fit model
options = optimset('fminsearch');  % specifying Simplex routine settings
options = optimset(options,'Display','Iter');  % display progress to command line

disp('Fitting model...'), tic
[pm] = fminsearchbnd(@(pm) find_optimal_bound(pm,sigmaM,sigmaD,lambda,theta,stimIn,delays,deltas),inits,lb,ub,options);  % parameters: [sigma bound]
toc

optimal_bounds(sets).noise = sigmaM; 
optimal_bounds(sets).optimal = pm; 

end 

save(['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/optimal_criterion/optimal_bounds.mat'] , 'optimal_bounds'); 