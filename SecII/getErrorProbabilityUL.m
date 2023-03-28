function [avg_error, reg_ul] = getErrorProbabilityUL(s, n, rho, rate, g_list, ghat_list, sigma_sq_list)
nbrOfRealizations = length(g_list);
eps_ul = nan(1, nbrOfRealizations);

for i = 1:nbrOfRealizations
    g = g_list(i);
    ghat = ghat_list(i);
    sigma_sq= sigma_sq_list(i);
    
    betaA_ul = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_ul = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_ul = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_ul*betaB_ul);
    
    
    preterm_ul = log(1+s*rho * abs(ghat)^2);
    if sum(isnan(betaA_ul))~= 0
        keyboard
    end
    [eps_ul(i), reg_ul(i) ] = saddlepoint_approximation(n, rate, betaA_ul, betaB_ul, nu_ul, preterm_ul);
    
end

avg_error=mean(eps_ul);

end