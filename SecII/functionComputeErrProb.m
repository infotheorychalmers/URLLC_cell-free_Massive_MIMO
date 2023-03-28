function [eps_ul,  s_val_ul] = functionComputeErrProb(H, n, rate, eps_target, combiner)

%power
pow = 1;

%Compute MR combining
V_MR = H;

if strcmp(combiner,'MR')
    
    % Combiner
    v1 = V_MR(:,1);
    
elseif strcmp(combiner,'MMSE')

    %Compute MMSE combining
    V_MMSE = ((V_MR*V_MR')+eye(size(H,1)))\V_MR;
    
    % Combiner
    v1 = V_MMSE(:,1);
    
end
v1 = v1./vecnorm(v1);

%Computation of g
g = v1' * H(:,1);

%Computation of noise deviation standard (including interference)
sigma_sq = v1'*v1 + sum(abs(v1'*H).^2) - abs(v1' * H(:,1))^2;

%-------------------------------------
% EVALUATE RCUs for perfect CSI
f_ul = @(s) getErrorProbabilityUL(s, n, pow, rate, g, g, sigma_sq);

[eps_ul,  s_val_ul] = golden_search(f_ul,1e-8,1,1e-8,eps_target);