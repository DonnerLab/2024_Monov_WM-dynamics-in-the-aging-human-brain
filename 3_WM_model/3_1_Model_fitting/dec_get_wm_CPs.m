% Recovers "different" response probabilities for a given combination of sigma
% and bound parameters (pm), for specified trials (data)
%
% Peter Murphy & Gina Monov, UKE Hamburg, 2022

function CPs = dec_get_wm_CPs(pm)
global modelspec
global model_log
global data_in
% Set simulation parameters

    sigmaM = []; %Initialize all possible free parameters
    sigmaD = [];
    bound = [];
    theta = [];
    lambda = [];
    % only feed in the ones specified in modelspec
    for r=1:length(model_log(model_log==1))
      
         if isempty(sigmaM) & ismember('sM',modelspec)
         sigmaM = pm(r);
         elseif isempty(sigmaD) & ismember('sD',modelspec)
         sigmaD = pm(r);
         elseif isempty(bound) & ismember('bound',modelspec)
         bound = pm(r);
         elseif isempty(theta) & ismember('theta',modelspec)
         theta = pm(r);
         elseif isempty(lambda) & ismember('lambda',modelspec)
         lambda = pm(r);
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


% Construct vector of model's "different" response probabilities to match the
% vector of participant responses
resp = data_in(:,3);  % participant responses

CPs = zeros(size(resp));
for d = 1:length(delays)
    for s = 1:length(deltas)
        CPs(data_in(:,1)==delays(d) & data_in(:,2)==deltas(s)) = cp(d,s);
    end
end