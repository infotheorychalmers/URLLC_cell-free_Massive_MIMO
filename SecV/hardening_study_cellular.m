function data = hardening_study_cellular(L, K, M_tot, infobits, n, np, rho_dBm, ESTIMATOR, COMBINER, nbrOfRealizations, antennaType, net_type, simNr, UE_CSI,CSI_AP,UE_POS)
vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix

DEBUG = 0;
rng('shuffle')

if DEBUG == 1
    %FLAGS
    
    K = 1; % Number of UEs in the network
    net_type = 'CELLULAR'; % type of network [CELLULAR,CENTRALIZED]. Only cellular here.
    L = 1; % Number of APs in the small-cell network
    M_tot = 144; % Number of antennas per AP
    np = 40; %number of pilots
    n = 300; %total number of channle uses
    infobits = 8*20; %transmission rate
    rho_dBm = -10; %transmit power [dBm]
    nbrOfRealizations=1e4; %number of saddlepoint realizations
    
    UE_POS = 6; % position of the UE (scaled by squareLength) - 1:(0.75 + 1i*0.5); 2:(0.775 + 1i*0.5); 3:(0.75 + 1i*0.75); 4:(0.775 + 1i*0.775); 5:(0.5 + 1i*0.75); 6:(0.5 + 1i*0.775)
    COMBINER = 'MR'; %what combiner to use [MR, M-MMSE, RZF]
    ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
    antennaType = 'ULA';
    
    simNr = 1;
    UE_CSI=0; % UE_CSI=1 means same CSI at the BS and UE
    CSI_AP = 0;
end
if UE_CSI ~= 0
    n_ul = (n - 2*np) / 2; %length of UL phase
    n_dl = n - 2*np - n_ul; %length of DL phase
else
    n_ul = (n - np) / 2; %length of UL phase
    n_dl = n - np - n_ul; %length of DL phase
end
M_list = (2:sqrt(M_tot)).^2;
eps_target = 1e-5;
squareLength = 150; % size of scenario (meter)
ASDdeg = 25; %Angular standard deviation around the nominal angle (measured in degrees)
p = 10^(rho_dBm/10);

%% Run Simulation
hardening = nan(1,length(M_list));

% Initializations:
g_ul = nan(nbrOfRealizations,1);
g_dl = g_ul;
ghat_dl = g_ul;
sigma_sq_dl = g_ul;

eps_dl = nan(1,length(M_list));
s_val_dl = nan(1,length(M_list));
tic
for mm = 1:length(M_list)
        M = M_list(mm);
        
        %% Generate setup and compute channel statistics:
        [R_AP,pilotIndex,~,UEpositions,APpositions] = generateSetup_hardening_study(L,K,M,np,ASDdeg,UE_POS,squareLength);
        %R_values = squeeze(R_AP) % Uncomment to display the correlation matrices
        k = randi(K,1); %pick a random UE
        % Create batches:
        batchsize = min(1e2, nbrOfRealizations);
        batches = nbrOfRealizations / batchsize;
        
        %Need to do this for memory issues
        for batch_indx =1:batches
            
            %% Generate channel and channel estimates
            [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R_AP,batchsize,L,K,M,np,pilotIndex,p,ESTIMATOR);
            % MSE = trace(C(:,:,1,1,1))/trace(Rscaled(:,:,1,1,1));
            %hardening_prev_combining (ll,mm) = var(vecnorm(H_AP).^2)/mean(vecnorm(H_AP).^2)^2;
            %% Generate combiner and precoding matrices for all UEs in each cell
            if CSI_AP == 0
                [~, W] =  functionCombinerPrecoder(Hhat_AP,B_AP,R_AP,batchsize,M,K,L,p, COMBINER, net_type);
            else
                [~, W] =  functionCombinerPrecoder(H_AP,B_AP,R_AP,batchsize,M,K,L,p, COMBINER, net_type);
            end
            for batch_sample = 1:batchsize
                
                combinedChannels = H_AP(:,batch_sample,k)' * squeeze(W(:,batch_sample,:));
                
                g_dl((batch_indx-1)*batchsize + batch_sample) = combinedChannels(k);
                if UE_CSI == 2
                    ghat_dl((batch_indx-1)*batchsize + batch_sample) = functionChannelEstimatesDL(combinedChannels,length(batch_sample),np,pilotIndex,p,ESTIMATOR,k);
                end
                combinedChannels(k) = []; %remove the UEs channel to only have noise terms left
                
                sigma_sq_dl((batch_indx-1)*batchsize + batch_sample) = 1 +  p*sum(abs(combinedChannels).^2);
            end
            if UE_CSI ~= 2
                ghat_dl((batch_indx-1)*batchsize + 1 : batch_indx*batchsize) = sum(conj(H_AP(:,:,k)) .* W(:,:,k), 1); % this will be used to estimate the channel hardening for DL
            end
            
            
            
        end
        %check channel hardening
        hardening(mm) = var(abs(ghat_dl).^2)/mean(abs(ghat_dl).^2)^2;
        %             g_test = ghat_dl;
        %             ghat_dl = mean(ghat_dl) * ones(size(g_dl)); %estimate channel hardening
        %             figure; plot(abs(g_test).^2,'b'); hold on;plot(1:length(g_test),abs(ghat_dl).^2,'r')
        
        time = toc;
        disp(['Scenario + estimation + combiner takes: ' num2str(time) ' seconds']);
        
        if UE_CSI == 0
            ghat_dl = mean(ghat_dl) * ones(size(g_dl)); %estimate channel hardening
        end
        
        %-----------------------------------------------
        % Estimate the uplink error probability
        %-----------------------------------------------
%         tic
%         
%         rate = infobits / n_ul;
%         f_ul = @(s) getErrorProbabilityUL(s, n_ul, p, rate, g_ul, ghat_ul, sigma_sq_ul);%only a tenth of samples used to find candidate
%         [eps_ul,  s_val_ul] = golden_search(f_ul,1e-8,1,1e-8,eps_target);
%         time2 = toc;
%         
%         disp(['UL search took: ' num2str(time2) ' seconds'])
        
        %-----------------------------------------------
        % Estimate the downlink error probability
        %-----------------------------------------------
        rate = infobits / n_dl;
        
        tic
        f_dl = @(s) getErrorProbabilityDL(s,n_dl, p, rate, g_dl, ghat_dl, sigma_sq_dl); %only a tenth of samples used to find candidate
        [eps_dl(mm),  s_val_dl(mm)] = golden_search(f_dl,1e-8,1,1e-8,eps_target);
        time4=toc;
        disp(['DL search took: ' num2str(time4) ' seconds'])
        
        %disp(['UL: epsilon = ' num2str(eps_ul(ll,mm)) ' and s = ' num2str(s_val_ul(ll,mm))]);
        disp(['DL: epsilon = ' num2str(eps_dl(mm)) ' and s = ' num2str(s_val_dl(mm))]);
        fprintf('\n')
        

    if eps_dl(mm) < eps_target/10
        break
    end
end

%-------------------------------------------------
%       FIGURES
%-------------------------------------------------

% Plot variance channel hardening after combining
figure
plot(M_list, hardening)
title(['cellular network with L = ' num2str(L) ', AP CSI = ' num2str(CSI_AP) ', UE CSI = ' num2str(UE_CSI)])
xlabel('total number of antennas in the system, LM')
ylabel('Normalized variance hardening metric')
grid on
% Plot downlink average error probability
figure
semilogy(M_list,eps_dl)
ylim([1e-5,1])
title(['cellular network with L = ' num2str(L) ', AP CSI = ' num2str(CSI_AP) ', UE CSI = ' num2str(UE_CSI)])
xlabel('total number of antennas in the system, LM')
ylabel('downlink average error probability')
grid on

time = toc;
disp(['Scenario + estimation + combiner takes: ' num2str(time) ' seconds']);


%-------------------------------------------------
%       SAVE FILE
%-------------------------------------------------

data.ASDdeg = ASDdeg;
data.L = L;
data.K=K;
data.M_list=M_list(M_list<=M);
data.np=np;
data.n = n;
data.rate =  infobits/n_ul;
data.rho_dBm = rho_dBm;
data.nbrOfRealizations = nbrOfRealizations;
data.COMBINER = COMBINER;
data.ESTIMATOR = ESTIMATOR;
data.hardening = hardening(~isnan(eps_dl));
data.s_val_dl = s_val_dl(~isnan(eps_dl));
data.avg_error_dl = eps_dl(~isnan(eps_dl));
data.CSI_AP = CSI_AP;
data.CSI_UE = UE_CSI;
data.UE_POS = UE_POS;
data.UEpositions = UEpositions;
data.APpositions = APpositions;
data.squareLength = squareLength;
% %-------------------------------------------------
% %       SAVE FILE
% %-------------------------------------------------

filename = ['hardening_simNr_' num2str(simNr) '_' net_type '_' ESTIMATOR '_' COMBINER '_' antennaType '_M_' num2str(M_tot) '_L_' num2str(L) '_CSI_AP_' num2str(CSI_AP) '_CSI_UE_' num2str(UE_CSI) '_K_' num2str(K) '_UE_POS_' num2str(UE_POS) '_rho_' num2str(rho_dBm) '_np_' num2str(np) '_n_' num2str(n) '.mat']

save(filename, 'data', '-v7.3')

disp('Finished')

end


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

function [avg_error, reg_dl] = getErrorProbabilityDL(s, n, rho, rate, g_list, ghat_list, sigma_sq_list)


nbrOfRealizations = length(g_list);
%check channel hardening
%  g_test = sum(conj(H(:,:,k,j,j)) .* squeeze(W(:,:,:,j)),1);
%  plot(abs(g_test),'b'); hold on;plot(1:length(g_test),ones(size(g_test))*abs(g_hat_dl),'r--')

eps_dl = nan(1, nbrOfRealizations);
for i = 1:nbrOfRealizations
    
    g = g_list(i);
    ghat = ghat_list(i);
    sigma_sq= sigma_sq_list(i);
    
    betaA_dl = s*rho*abs(g-ghat)^2 + s*sigma_sq;
    betaB_dl = s*(rho*abs(g)^2 + sigma_sq) / (1+s*rho*abs(ghat)^2);
    
    sigma_v = abs(g)^2 *rho + sigma_sq;
    gamma = s/(1 + s*rho*abs(ghat)^2);
    nu_dl = s*gamma*abs(sigma_v - rho* g'*ghat)^2 / (betaA_dl*betaB_dl);
    
    
    preterm_dl = log(1+s*rho * abs(ghat)^2);
    %----------------------------
    
    
    [eps_dl(i), reg_dl(i)] = saddlepoint_approximation(n, rate, betaA_dl, betaB_dl, nu_dl, preterm_dl);
    
end

avg_error=  mean(eps_dl);

end

