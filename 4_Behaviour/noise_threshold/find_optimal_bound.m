% Find optimal bound at a given amount of memory noise 
% Gina Monov, UKE, 2023

function SSR = find_optimal_bound(pm,sigmaM,sigmaD,lambda,theta,stimIn,delays,deltas)

% Set simulation parameters
bound = pm(1);

% Generate single-trial choice probabilities and choices
for t = 1:size(stimIn,1)
    cp(t,1) = run_wd_DecRule(sigmaM,sigmaD,bound,theta,delays(t),deltas(t));  % CP given memory noise and decision rule
    if ~isempty(lambda)
       cp(t,1) = cp(t).*(1-(1-exp(-lambda*delays(t)))) + 0.5.*(1-exp(-lambda*delays(t)));  % add memory lapse component if desired
    end
end

cp(cp>1) = 1; % adjusting CPs that slightly exceed 1 (occurs rarely due to numerical imprecisons) 

resp(:,1) = cp; 

% Extract correct (1) and incorrect (0)

for tt = 1:length(resp)
    if deltas(tt) == 0
      resp(tt) = 1-resp(tt); 
    end 
end 

SSR = 1-mean(resp); 