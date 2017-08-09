function [ output_args ] = mu_shift( input_args )
%BIG SUM
%   Lil details.

%Find indicies where the difference between MU for a selected phase-arc is
%greater than the MU difference threshold.
mu_diff_indicies = find(diff( mu_e(sorted_phase_arc{1,1}) ) > mu_diff_threshold);

for i = 1:length(mu_diff_indicies)
    mu_initial = sorted_phase_arc{1,1}(mu_diff_indicies(i));
    


end

