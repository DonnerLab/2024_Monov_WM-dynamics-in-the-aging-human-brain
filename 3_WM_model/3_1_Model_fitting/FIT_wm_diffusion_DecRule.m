% Fits Wiener diffusion model with sigmoidal decision rule to behaviour
% on delayed-match-to-sample working memory task. Five free parameters:
% noise in the diffusion process (sigmaM); decision noise deterimning slope
% of sigmoidal decision rule (sigmaD); the criterion for "non-match" 
% judgements (denoted 'bound', determining inflection point of decision
% function) applied to the comparison of the memory trace and target
% spatial positions; the percentage of static lapses (theta) capturing
% sensory encoding or motor exectution lapses; and the hazard rate of a
% memory lapse (lambda), capturing hazard rate with which memory trace is
% lost (fixed lambda = growing hazard function over time). All forms of
% lapse produce random responses. Model is fit through Particle Swarm Optimization
% Birge B (2003) PSOt - a particle swarm optimization toolbox for use with Matlab. 
% In: Proceedings of the 2003 IEEE Swarm Intelligence Symposium. SIS?03 (Cat. No.03EX706), pp 182?186. 
% Indianapolis, IN, USA: IEEE. Available at: http://ieeexplore.ieee.org/document/1202265/ [Accessed November 7, 2022].
%
% Peter Murphy, Trinity College Dublin & Gina Monov, UKE, 2022

function[pm] = FIT_wm_diffusion_DecRule(n)
mode = 1; % set mode either to 1 (MCI study) or 2 (YHC) or 3 (syntehtic datasets for parameter recovery)

% Selct which set of free parameters to fit 
global modelspec % modelspec defines which parameters will be included 
global model_log
%modelspec = 'sM+bound'; 
%modelspec = 'sM+bound+theta'; 
%modelspec = 'sM+bound+lambda'; 
%modelspec = 'sM+bound+theta+lambda'; 
%modelspec = 'sM+sD+bound'; 
modelspec = 'sM+sD+bound+theta+lambda'; 
model_log = zeros(1,5); % logical of free parameters in specified model
  
% load details
if mode == 1
    loadpath = '/mnt/homes/home028/gmonov/behav_eye_data/';
    load('/mnt/homes/home028/gmonov/meg_analysis/datasets_overview.mat')
    allsubj={};
    for f = 1:length(datasets_overview)
        if datasets_overview(f).sess == 1
            allsubj = horzcat(allsubj,datasets_overview(f).subj);
        end
    end
    
elseif mode == 2
    loadpath = '/mnt/homes/home028/gmonov/SCZ/Data/working_memory/';
    load(['/mnt/homes/home028/gmonov/SCZ/Subjects/scz_wm_allsubj.mat']);
    allsubj = scz_wm_allsubj;
    
elseif mode == 3
    loadpath = ['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/gen_datasets/',modelspec, '/']; 
    allsubj={};
    for z = 1:100
        allsubj = horzcat(allsubj,num2str(z)); % create cell with dataset names 1 to 100 
    end 
end

% Set global variables to pass to PSO algorithm
global data_in

subj = allsubj{n};


addpath(genpath('/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/particle_swarm_WM'))


% Define parameter estimation settings
if ismember('sM',modelspec)
range_p.sigmaM = [0.1 30];  % parameter bounds
model_log(1) = 1; 
end 
if ismember('sD',modelspec)
range_p.sigmaD = [1e-7 30];
model_log(2) = 1; 
end
if ismember('bound',modelspec)
range_p.bound  = [0 40];
model_log(3) = 1; 
end
if ismember('theta',modelspec)
range_p.theta  = [0 0.49];
model_log(4) = 1; 
end
if ismember('lambda',modelspec)
range_p.lambda = [0 0.8];
model_log(5) = 1; 
end

mv = [0.2;...  
    % maximum particle velocities
      0.1;...
      0.2;...
      0.001;...
      0.001]';
  
mv = mv(model_log==1); % remove free parameters that are not included  

seeds.sigmaM = [5 1];  % good seed distributions for a selection of particles - [mean sd]
seeds.sigmaD = [0.5 0.4];
seeds.bound  = [10 3];
seeds.theta  = [0.02 0.04];
seeds.lambda = [0.02 0.04];

% Seed random number generator
seed = round(sum(100*clock)); %never the same seed
rand('state', seed);


if mode ~= 3 % only pull behavior and make exceptions for real data 
% Load behavioural data and concatenate
all_behav = [];
if mode == 1
    fullpath = [loadpath,filesep,subj,filesep,'S1',filesep,'Behaviour'] % for MCI
elseif mode ==2
    fullpath = [loadpath,subj,'/S2/Behaviour'] %for young subjects
end
files = dir([fullpath,filesep,'*.mat']);

% Excluding 1st block
if strcmp(subj,'01')||strcmp(subj,'25') ||strcmp(subj,'02') % terminated after a few (2 and 9) trials / block ecxlusion due to performance 
    files(1)=[];
end


for f = 1:length(files)
    if mode == 1
        load([loadpath,filesep,subj,filesep, 'S1', filesep,'Behaviour',filesep, files(f).name]) %MCI
    elseif mode == 2
        load([fullpath,filesep,files(f).name]) % YHC
    end
    all_behav = [all_behav; Behav];
end

% Excluding trials with too long RT/ too short RT/ wrong button press
R = all_behav(:,7);
all_behav(R>=(mean(R)+4*std(R)),:) = [];
R = all_behav(:,7);
all_behav(R<=0.2,:) =[];
choice = all_behav(:,5);
all_behav(choice==99,:)=[];

% Pull delay durations, memorandum-target location differences, and response
delays = all_behav(:,3);
deltas = round(abs(all_behav(:,2)-all_behav(:,1)),1);  % rounding because otherwise there are small rounding errors which make indexing annoying
resp = all_behav(:,5);
resp = resp.*-1;  % flipping sign of resp vector so that 1="non-match" response, -1="match" response
resp(resp==-1) = 0; % 0 = "match", 1 = "non-match" (this accords w/ decision rule, which determines p("non-match"))


else 
    load([loadpath,subj,'.mat']) % load generated data set 
    delays = dat(:,1);
    deltas = dat(:,2); 
    resp = dat(:,3); 
    
end 
% Defining PSO options (see pso_Trelea_vectorized.m for details)
P(1)=0;  P(2)=1000;     P(3)=100;    P(4:13)=[1.6 1.9 0.9 0.4 400 1e-25 250 NaN 0 1];
% display  n_iterations  n_particles       acceleration, inertia, tolerance, etc

% Seeding first n particles with parameters drawn from realistic distributions
n_seeded = 50;
PSOseedValue=[];
if ismember('sM',modelspec)
PSOseedValue(1:n_seeded,end+1) = seeds.sigmaM(1)+(randn(n_seeded,1).*seeds.sigmaM(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.sigmaM(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.sigmaM(1)),end) = range_p.sigmaM(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.sigmaM(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.sigmaM(2)),end) = range_p.sigmaM(2); end
end

if ismember('sD',modelspec)
PSOseedValue(1:n_seeded,end+1) = seeds.sigmaD(1)+(randn(n_seeded,1).*seeds.sigmaD(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.sigmaD(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.sigmaD(1)),end) = range_p.sigmaD(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.sigmaD(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.sigmaD(2)),end) = range_p.sigmaD(2); end
end

if ismember('bound',modelspec)
PSOseedValue(1:n_seeded,end+1) = seeds.bound(1)+(randn(n_seeded,1).*seeds.bound(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.bound(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.bound(1)),end) = range_p.bound(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.bound(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.bound(2)),end) = range_p.bound(2); end
end 

if ismember('theta',modelspec)
PSOseedValue(1:n_seeded,end+1) = seeds.theta(1)+(randn(n_seeded,1).*seeds.theta(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.theta(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.theta(1)),end) = range_p.theta(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.theta(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.theta(2)),end) = range_p.theta(2); end
end 

if ismember('lambda',modelspec)
PSOseedValue(1:n_seeded,end+1) = seeds.lambda(1)+(randn(n_seeded,1).*seeds.lambda(2));
if ~isempty(find(PSOseedValue(:,end)<range_p.lambda(1))), PSOseedValue(find(PSOseedValue(:,end)<range_p.lambda(1)),end) = range_p.lambda(1); end
if ~isempty(find(PSOseedValue(:,end)>range_p.lambda(2))), PSOseedValue(find(PSOseedValue(:,end)>range_p.lambda(2)),end) = range_p.lambda(2); end
end 

% Concatenating parameter ranges
par_range = [];
if ismember('sM',modelspec)
   par_range = [par_range;range_p.sigmaM];
end
if ismember('sD',modelspec)
   par_range = [par_range;range_p.sigmaD];
end
if ismember('bound',modelspec)
   par_range = [par_range;range_p.bound];
end
if ismember('theta',modelspec)
   par_range = [par_range;range_p.theta];
end
if ismember('lambda',modelspec)
   par_range = [par_range;range_p.lambda];
end


% Randomly seeding remaining particles within prespecified bounds
PSOseedValue(size(PSOseedValue,1)+1:P(3),1:size(PSOseedValue,2)) = normmat(rand([P(3)-n_seeded,size(PSOseedValue,2)]),par_range',1);

% Running PSO routine
data_in = [delays deltas resp];
disp('Fitting model...'), tic
[output,te,tr] = pso_Trelea_vectorized_WM('wm_diffusion_DecRule_CE',length(mv),mv,par_range,0,P,'goplotpso',PSOseedValue);
toc

% Store parameter estimates
pm = output(1:end-1);
ce = output(end);


% Retrieve model's response probabilities for visualization of fit

CPs = dec_get_wm_CPs(pm);
CPs(CPs > 1 & CPs < 1.0001) = 1; % adjusting CPs that slightly exceed 1 (occurs rarely due to numerical imprecisons) 
if mode == 2
    save(['/mnt/homes/home028/gmonov/SCZ/WM_modeling/Dec_Rule/',modelspec,filesep,subj,'.mat'] , 'CPs','pm','all_behav','ce','modelspec');

elseif mode == 1
    save(['/mnt/homes/home028/gmonov/Modelling/Dec_Rule/',modelspec,filesep,subj,'.mat'] , 'CPs','pm','all_behav','ce','modelspec');
    
elseif mode == 3
    save(['/mnt/homes/home028/gmonov/Modelling/Dec_Rule_Modelling/gen_datasets/',modelspec,filesep,'Fits',filesep,subj,'.mat'] , 'CPs','pm','dat','ce','modelspec');
end