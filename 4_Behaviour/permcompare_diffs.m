% Compare difference of a difference using two-tailed permutation test 
% Gina Monov, UKE, Hamburg, 2023 

function comp_diffs_results=permcompare_diffs(a1,a2,b1,b2,nperms) 
% Specify inputs for group a and b (variable 1, variable 2) and number of
% permutations 
% function tests significance of the difference of group differences
% between two variables 

comp_diffs_results = table; 

% Data normalization
n_a = length(a1); 
n_b = length(b1); 

all_1 = [a1;b1]; 
all_2 = [a2;b2]; 

all_1 = zscore(all_1); 
all_2 = zscore(all_2); 

clear a1 a2 b1 b2 

a1 = all_1(1:n_a);
b1 = all_1(n_a+1:n_a+n_b);

a2 = all_2(1:n_a);
b2 = all_2(n_a+1:n_a+n_b);

comp_diffs_results.group_difference_1 = mean(a1)-mean(b1); %Compute group difference for variable 1
comp_diffs_results.group_difference_2 = mean(a2)-mean(b2); %Compute group difference for variable 2

comp_diffs_results.delta_difference = comp_diffs_results.group_difference_1-comp_diffs_results.group_difference_2; % delta of the group differences 


% Concatenate both groups 

merged_data_1 = [a1;b1]; 
merged_data_2 = [a2;b2]; 

    for l = 1:nperms
        sort_ids = zeros(n_a+n_b,1); 
        sort_ids(1:n_a,1) = 1; 
        idx = randperm(n_a+n_b); 
        sort_ids = sort_ids(idx); % shuffle group allocation
        test_group_difference_1 = mean(merged_data_1(sort_ids==1))-mean(merged_data_1(sort_ids==0)); 
        test_group_difference_2 = mean(merged_data_2(sort_ids==1))-mean(merged_data_2(sort_ids==0)); 
        
        null_rs(l) = test_group_difference_1-test_group_difference_2; 
        clear test_*
    end 
   
  % Compute percentile to obtain p-value for the null-hypothesis 
    nulls = null_rs(1:end);  
    if comp_diffs_results.delta_difference>0
       p_v = length(nulls(nulls > comp_diffs_results.delta_difference))./nperms.*2; % two-tailed
    elseif comp_diffs_results.delta_difference<0
       p_v = length(nulls(nulls < comp_diffs_results.delta_difference))./nperms.*2; % two-tailed
    end 
   
   comp_diffs_results.diff_p = p_v; 
   comp_diffs_results.nulls{1} = nulls; 

end 