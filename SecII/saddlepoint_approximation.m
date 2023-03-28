function [eps_out, region] = saddlepoint_approximation(n, R, a, b, c, d)
%n is the blocklength
%R is in bits
%a = betaA;
%b = betaB;
%c = nu;
%d = preterm;

R=R*log(2); %rate in nats

% Useful definitions for compactness:
xi      = 2*a*b*(1-c);
lambda  = sqrt((a-b)^2+4*a*b*(1-c));
gamma   = a - b;

%Region of convergence:
zeta_low = -(gamma + lambda)/xi;
zeta_high = (lambda - gamma)/xi;


% CGF and derivatives as defined in the paper:
kappa0 = @(t) -t*d - log(1 + (b-a)*t - a*b*(1-c)*t.^2);
kappa1 = @(t) -(-2*a*b*(1-c)*t +b-a)./(-a*b*(1-c).*t.^2+t*(b-a)+1) - d;
kappa2 = @(t) (b-a-2*a*b*(1-c)*t).^2./(1+(b-a)*t-a*b*(1-c)*t.^2).^2 + (2*a*b*(1-c))./(1+(b-a)*t-a*b*(1-c).*t.^2);

% Some useful functions:
Rcr = -kappa1(1); %critical rate
Is = -kappa1(0); %generalized MI


Psi  = @(u,k2) exp(n*u^2*k2/2 + log(qfunc(u*sqrt(n*k2))));
PsiTilde = @(x,y)  exp(n*x*(Rcr-R+kappa2(1)/2) + log(qfunc(x*sqrt(n*kappa2(1)) + y*n*(Rcr - R)/sqrt(n*kappa2(1)))));

% %--------------------------------------------
% %find start value for stationary equation
% small_tmp = eps;
% zeta_candidate_low = zeta_low+small_tmp;
% while kappa1(zeta_candidate_low) > kappa1(0) || abs(kappa1(zeta_candidate_low))== inf%because of numerical issues.
%     small_tmp= 2*small_tmp;
%     zeta_candidate_low = zeta_low+small_tmp;
% end
% while -kappa1(zeta_candidate_low) > R && zeta_candidate_low < zeta_high%search while we are still above R
%     zeta_candidate_low = zeta_candidate_low + 1;
% end
% zeta_candidate_low = max(zeta_candidate_low-2, zeta_low); 
% 
% %find end and start value for stationary equation -kappa'(zeta)=R
% small_tmp = eps;
% zeta_candidate_high = zeta_high-small_tmp;
% while kappa1(zeta_candidate_high) == inf ||  kappa1(zeta_candidate_high) < kappa1(0) %because of numerical issues.
%     small_tmp= 2*small_tmp;
%     zeta_candidate_high = zeta_high-small_tmp ;
% end
% while -kappa1(zeta_candidate_high) < R && zeta_candidate_high > zeta_low%search while we are still below R
%     zeta_candidate_high = zeta_candidate_high - 1;
% end
% zeta_candidate_high = min(zeta_candidate_high+2, zeta_high);
% 
% if abs(zeta_candidate_high -zeta_high) < 1 ||  abs(zeta_candidate_low - zeta_low) < 1 % we need to look very close to the ROC border
%     interval = linspace(zeta_candidate_low, zeta_candidate_high, 1e4) ;
% else
%     interval =linspace(zeta_candidate_low, zeta_candidate_high, 1e3) ;
% end
% 
% %--------------------------------------------
% 
% 
% 
% zeta_candidates = interval;
% zeta_candidates(1)=[];
% zeta_candidates(end)=[];
% k1 = kappa1(zeta_candidates);
% 
% %make sure we have no infs because we got too close to RoC endpoints
% [~,iduniq] = unique(k1);
% k1 = k1(iduniq);
% zeta_candidates= zeta_candidates(iduniq);
% idinf = find(isinf(k1));
% k1(idinf) = [];
% zeta_candidates(idinf)= [];
% 
% if sum(isinf(k1)) ~= 0 || sum(isnan(k1)) ~= 0 || sum(isinf(zeta_candidates)) ~= 0 || sum(isnan(zeta_candidates)) ~= 0
%    keyboard 
% end
% zeta = interp1(k1, zeta_candidates, -R); % Finding zeta s.t. kappa'(zeta) = threshold.

% Find the zeta satisfying the stationary equation:
func_kappa1_R = @(t) kappa1(t) + R;
error_zeta = 1;
iter = 0;
small_tmp = eps;
while error_zeta>1e-6
    if iter>20
        if sign(func_kappa1_R(zeta_interval(1)))==-1 &&  sign(func_kappa1_R(zeta_interval(2))) ==-1
            disp('Numerical issues! I cant find a zeta satisfying the stationary equation: R>Is for every zeta candidate. I will set epsilon = 1.')
            eps_out =1;
            region=1;
        elseif sign(func_kappa1_R(zeta_interval(1)))==1 &&  sign(func_kappa1_R(zeta_interval(2))) ==1
            disp('Numerical issues! I cant find a zeta satisfying the stationary equation: R< for every zeta candidate. I will set epsilon = 0.')
            eps_out =0;
            region=3;
        elseif isinf(func_kappa1_R(zeta_interval(1))) || isinf(func_kappa1_R(zeta_interval(2)))
            disp('Numerical issues! Something weird has happened when exploring the stationary equation. I cant solve numerically the problem and I set epsilon =1.')
            eps_out=1;
            region=1;
        end
        break;
    end
    small_tmp=small_tmp*10;
    iter = iter+1;
    zeta_interval = [zeta_low+small_tmp, zeta_high-small_tmp]; %[zeta_low zeta_high];
    if isinf(func_kappa1_R(zeta_interval(1))) || isinf(func_kappa1_R(zeta_interval(2)))
        continue
    elseif sign(func_kappa1_R(zeta_interval(1)))==-1 &&  sign(func_kappa1_R(zeta_interval(2))) ==-1
        continue
    elseif sign(func_kappa1_R(zeta_interval(1)))==1 &&  sign(func_kappa1_R(zeta_interval(2))) ==1
        continue
    end
    %options = optimset('Display','iter'); % show iterations
    [zeta, error_zeta] = fzero(func_kappa1_R,zeta_interval);
    error_zeta = abs(error_zeta);
end

if zeta > 1 % R < Rcr
    k0 = kappa0(1);
    eps_out = exp(n*(k0+R)) *(PsiTilde(1,1) + PsiTilde(0,-1));
    region = 3;
elseif zeta < 0 %R > Is
    k0 = kappa0(zeta);
    k2 = kappa2(zeta);
    error_exp   = exp(n*(k0+zeta*R)); % Error exponent.
    
    eps_out = 1 - error_exp * (Psi(-zeta,k2) - Psi(1-zeta, k2));
    
    region = 1;
else %Rcr < R < Is
    k0 = kappa0(zeta);
    k2 = kappa2(zeta); 
    error_exp   = exp(n*(k0 + zeta*R)); % Error exponent.
    
    eps_out = error_exp * (Psi(zeta,k2) + Psi(1-zeta, k2));
    
    region = 2;
end
if isnan(eps_out)
    keyboard;
end
if ~isnan(eps_out)
    eps_out = min(1, eps_out);
end

end
