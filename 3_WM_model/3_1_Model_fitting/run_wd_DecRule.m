% Computes p("non-match") choice probabilities given Wiener diffusion
% process with specified noise term ('sigmaM') and subject to sigmoidal
% decision rule with inflection point 'bound', slope 'sigmaD' and symmetric
% asymptotes 'theta' and 1-'theta'. Decision variable is defined as
% proximity of memory representation to location of test stimulus.
% Individual choice probabilities are computed for different delay
% durations/memorandum-target differences. 
% 
% Peter Murphy, Trinity College Dublin & Gina Monov, UKE, 2022

function [p] = run_wd_DecRule(sigmaM,sigmaD,bound,theta,delays,deltas)

%if model is fitted without decision noise set scale parameter to zero 
if isempty(sigmaD) 
    sigmaD = 0; 
end 
%if model is fitted without theta set lapse parameter to zero 
if isempty(theta) 
   theta = 0; 
end 
% Values at which to compute memory representation PDFs
dx = 0.05;   % resolution of PDF
x_in = -360:dx:360; 

 
% Compute SD of memory representation distributions @ each probe time
sds = sqrt(delays.*(sigmaM.^2));

% Loop through delay durations/mem-target deltas and calculate "same" response probabilities
for d = 1:length(delays)
    for s = 1:length(deltas)
        % compute analytical PDF of memory representation @ this delay

         ndf = 1./(sqrt(2.*pi.*(sds(d).^2))).*dx; 
         denom = 2.*(sds(d).^2); 
         dvpdf = ndf.* (exp(-(((x_in(x_in>=0)-deltas(s)).^2)./denom))) + ndf.* (exp(-(((x_in(x_in>=0)+deltas(s)).^2)./denom))); 
         dvpdf = dvpdf./sum(dvpdf); %renormalizing to deal with edge cases 

        % compute choice probability integrating over full DV PDF
         p(d,s) = sum((theta + (1-2*theta)./(1+exp(-(x_in(x_in>=0)-bound)./sigmaD))).*dvpdf);
         clear dvpdf
    end
end