% Computes cross entropy between participant responses and model's response
% probabilities, the latter computed via simulated Wiener diffusion
% processes subject to sigmoidal decision rule.
% 
% Inputs: pm = parameters [sigmaM sigmaD bound theta lambda]
%         data = trials*data-type matrix (data-type: [delay distance resp])
% 
% Peter Murphy, Trinity College Dublin & Gina Monov, UKE, 2022

function CE = wm_diffusion_DecRule_CE(pm)

% Retrieving global variables from initializing script
global data_in
global modelspec
global model_log
% Looping through PSO particles
for p = 1:size(pm,1)
    
    % Set simulation parameters
    sigmaM = []; %Initialize all possible free parameters
    sigmaD = [];
    bound = [];
    theta = [];
    lambda = [];
    %only feed in the ones specified in modelspec
    for r=1:length(model_log(model_log==1))
      
         if isempty(sigmaM) & ismember('sM',modelspec)
         sigmaM = pm(p,r);
         elseif isempty(sigmaD) & ismember('sD',modelspec)
         sigmaD = pm(p,r);
         elseif isempty(bound) & ismember('bound',modelspec)
         bound = pm(p,r);
         elseif isempty(theta) & ismember('theta',modelspec)
         theta = pm(p,r);
         elseif isempty(lambda) & ismember('lambda',modelspec)
         lambda = pm(p,r);
        end 
     
    end 
    
    
    delays = unique(data_in(:,1));
    deltas = unique(data_in(:,2));
    
    % Compute choice probabilities for each trial type
    cp = run_wd_DecRule(sigmaM,sigmaD,bound,theta,delays,deltas);  % returns p("non-match") for each delay/delta
    
    % Add memory lapses occuring @ hazard rate lambda: p(mem_lapse by t) = 1-exp(-lambda*t)
    if ~isempty(lambda) % only perform this if lambda is a free parameter
        for d = 1:length(delays)
           cp(d,:) = cp(d,:).*(1-(1-exp(-lambda*delays(d)))) + (ones(size(cp(d,:))).*0.5).*(1-exp(-lambda*delays(d)));
        end
    end 
    
    % Construct vector of model's "non-match" response probabilities to match
    % the vector of participant responses
    resp = data_in(:,3);  % participant responses
    
    CPs = zeros(size(resp));
    for d = 1:length(delays)
        for s = 1:length(deltas)
            CPs(data_in(:,1)==delays(d) & data_in(:,2)==deltas(s)) = cp(d,s);
        end
    end
    CPs(CPs > 1 & CPs < 1.0001) = 1; %adjusting CPs that slightly exceed 1 (occurs rarely due to numerical imprecisons) 
    % Compute cross entropy between participants' responses and model's
    % response probabilities
    CPs(CPs==1) = 0.99999; CPs(CPs==0) = 0.00001;  % adjusting extreme CPs that lead to inf cross-entropy values
    CE(p,1) = -sum((ones(size(resp))-resp).*log(1-CPs)+(resp.*log(CPs)));
end