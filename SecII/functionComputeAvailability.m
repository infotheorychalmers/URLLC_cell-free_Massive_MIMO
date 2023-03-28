function availability = functionComputeAvailability(eps, eps_target)

% Compute network availability 
[cdf, x] = ecdf(eps);

% [~ , target_indx] = min(abs(x - eps_target));
target_indx = find(abs(x - eps_target)==min(abs(x - eps_target)),1,'last');
if max(x) < eps_target
    target_indx = length(cdf);
elseif min(x) > eps_target
    target_indx = 1;
end

availability = cdf(target_indx);

end