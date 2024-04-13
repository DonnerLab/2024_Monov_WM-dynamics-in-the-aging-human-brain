% Compare correlation coefficients of two groups using two-tailed
% permutation test
% Gina Monov, UKE, Hamburg, 2023 

function comp_corr_results=permcompare_corrs(a1,a2,b1,b2,nperms) 
% Specify inputs for group a and b (variable 1, variable 2) and number of
% permutations 
% OUTPUT: table consisting of correlation coeffiecients and corresponding
% p-values for group a and group b, delta-r and p-value for delta-r, and
% distribution of delta-r under the null hypothesis 
comp_corr_results = table; 

rows = 'Rows'; 
pw = 'Pairwise'; % can theoretically deal with nans
tt = 'type'; 
corr_type = 'Pearson'; % change correlation type if desired 

% Find real correlation coefficients and differences
[r_a,pval_a] = corr(a1,a2,rows,pw,tt,corr_type); 
[r_b,pval_b] = corr(b1,b2,rows,pw,tt,corr_type); 

comp_corr_results.Pearons_r_group_a = r_a;
comp_corr_results.p_group_a = pval_a;
comp_corr_results.Pearons_r_group_b = r_b;
comp_corr_results.p_group_b = pval_b;
delta_r = r_a-r_b; 

comp_corr_results.Delta_r = delta_r; 

% Concatenate both groups 
all_1 = [a1;b1]; 
all_2 = [a2;b2]; 

n_a = length(a1); 
n_b = length(b1); 

    for l = 1:nperms
        sort_ids = zeros(n_a+n_b,1); 
        sort_ids(1:n_a,1) = 1; 
        idx = randperm(n_a+n_b); 
        sort_ids = sort_ids(idx); % shuffle group allocation
        eval(['[r1,pval] = corr(all_1(sort_ids == 1),all_2(sort_ids == 1),rows,pw,tt,corr_type);']); % random group a
        eval(['[r2,pval] = corr(all_1(sort_ids == 0),all_2(sort_ids == 0),rows,pw,tt,corr_type);']); % random group b
        null_rs(l) = r1-r2; 
        clear r1 r2
    end 
   
    % Compute percentile to obtain p-value for the null-hypothesis 
    % (group correlations do not differ)
    nulls = null_rs(1:end);  
    if comp_corr_results.Delta_r>0
       p_v = length(nulls(nulls > comp_corr_results.Delta_r))./nperms.*2; % two-tailed
    elseif comp_corr_results.Delta_r<0
       p_v = length(nulls(nulls < comp_corr_results.Delta_r))./nperms.*2; % two-tailed
    end 
   
    comp_corr_results.Delta_p = p_v; 
    comp_corr_results.nulls{1} = nulls; 

end 