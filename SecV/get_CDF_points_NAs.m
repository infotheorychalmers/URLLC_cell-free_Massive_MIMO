function data = get_CDF_points_NAs(L, K, M, infobits, n, np, rho_dBm, ESTIMATOR, COMBINER, nbrOfRealizations,nbrOfPositions, antennaType, net_type, simNr, UE_CSI)
% Function get_CDF_points(L, K, M, infobits, n, np, rho_dBm, ESTIMATOR, COMBINER, nbrOfRealizations,nbrOfPositions, antennaType, net_type, simNr, UE_CSI): 
% Computes the average error probability in both UL and DL for randomly deployed scenarios.
% 
% The average error probabilities are estimated using the saddlepoint
% approximation (function saddlepoint_approximation.m). To do so, the
% function getErrorProbabilityUL(DL) computes the saddlepoint approximation
% for a given s. The parameter s is optimized using the function
% golden_search.m
%
% The scenario is generated using the function generateSetup.m
% The channel estimates are obtained using the function functionChannelEstimates.m
% The combiners/precoders are obtained using the function functionCombinerPrecoder.m
%
% This function contains three different network settings:
% CENTRALIZED: cell-free setting where the signal processing is performed at a CPU.
% DISTRIBUTED: cell-free setting where the channel estimation and spatial processing 
% is performed at the AP, and the CPU performs extra signal processing only
% based on the statistics of the channel and the results obtained at each
% AP after performing the spatial processing.
% CELLULAR: Classical Massive MIMO cellular network. It also includes the
% case of small cellular networks.
%
%
% INPUTS:
% L = number of cells
% K = number of users per cell
% M = number of antennas at BS
% infobits = information bits
% n = blocklength
% np = length pilot sequence
% rho_db = transmit power [dB]
% ESTIMATOR = Channel estimator to use [LS, MMSE]
% COMBINER = Combiner to use [MR, M-MMSE, RZF]
% nbrOfRealizations = number of realization in the Monte-Carlo simulations
% nbrOfPositions = Number of random positions to generate the CDF
% antennaType = ['ULA','UCA'] uniform linear array or uniform circular array
% simNO: Simulation number (optional). Used externally to generate more batches and get more points of the CDF
% UE_CSI => UE_CSI=0 means decoding relying on channel hardening
%           UE_CSI=1 means same CSI at the BS and UE
%           UE_CSI=2 means pilots in the DL


vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix

DEBUG = 0; % Set it to 1 to run a quick simulation using the parameters defined in the if condition below. 
rng('shuffle')
if nargin < 15
    UE_CSI=0;
end

if DEBUG == 1
    %FLAGS
    K = 40; % Number of UEs in the network
    L = 25; % Number of APs
    M = 16; % number of antennas per AP
    np = K; %number of pilots
    n = 300; %total number of channle uses
    infobits = 8*20; %transmission rate
    rho_dBm = -10; %transmit power [dBm]
    nbrOfRealizations=1e3; %number of saddlepoint realizations
    nbrOfPositions = 1;
    
    COMBINER = 'M-MMSE'; %what combiner to use [MR, M-MMSE, RZF]
    ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
    antennaType = 'ULA';
    net_type = 'DISTRIBUTED'; % type of network [CELLULAR,DISTRIBUTED,CENTRALIZED]
    
    simNr = 1;
    UE_CSI=2; % UE_CSI=0 means decoding relying on channel hardening
              % UE_CSI=1 means same CSI at the BS and UE
              % UE_CSI=2 means pilots in the DL
end
eps_target = 1e-5;

ASDdeg = 25; %Angular standard deviation around the nominal angle (measured in degrees)

if UE_CSI ~= 0
    n_ul = (n - 2*np) / 2; %length of UL phase
    n_dl = n - 2*np - n_ul; %length of DL phase
    
else
    n_ul = n/2 - np; %length of UL phase
    n_dl = n - n_ul - np; %length of DL phase
end
p = 10^(rho_dBm/10);
%% Run Simulation

% Initializations:
eps_ul = inf(nbrOfPositions,1);
eps_dl = inf(nbrOfPositions,1);
s_val_ul = inf(nbrOfPositions,1);
s_val_dl = inf(nbrOfPositions,1);

% Compute the perfromance for nbrOfPositions scenarios (randomly deployed users within the scenario)
parfor pos_idx = 1:nbrOfPositions % Go through all setups
    disp([num2str(pos_idx) ' setups out of ' num2str(nbrOfPositions)]);%Output simulation progress
    % Initializations:
    g_ul = nan(nbrOfRealizations,1);
    g_ul_all_dist = nan(L,nbrOfRealizations,K);
    ghat_ul = g_ul;
    sigma_sq_ul = g_ul;
    g_dl = g_ul;
    ghat_dl = g_ul;
    sigma_sq_dl = g_ul;
    v_k_samp = []; % Needed to avoid warning due to parfor
    
    tic
    %% Generate setup and compute channel statistics:
    [R_AP,pilotIndex,gainOverNoisedB_AP] = generateSetup(L,K,M,np,ASDdeg);
    
    %%%%%%%%%%%%%%%%%%%%%%% Only used in CELLULAR case %%%%%%%%%%%%%%%%%%%%
    % We create the masks to implement Level 1 of cooperation (small CELLULAR network where APs work as independent BSs)
    APassignment = zeros(M*size(gainOverNoisedB_AP,1),size(gainOverNoisedB_AP,2));
    [~,servingAP_cellular] = max(gainOverNoisedB_AP,[],1); % The AP with larger channel gain over noise serves the given UE.
    for ii =1:K
        APassignment(M*(servingAP_cellular(ii)-1)+1:M*(servingAP_cellular(ii)),ii) = 1; % Select the N antennas from the serving AP.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    k = randi(K,1); %pick a random UE
    
    % Create batches (We need to do this for memory issues):
    batchsize = min(1e2, nbrOfRealizations);
    batches = nbrOfRealizations / batchsize;
    
    for batch_indx =1:batches
        %% Generate channel and channel estimates
        [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R_AP,batchsize,L,K,M,np,pilotIndex,p,ESTIMATOR);
        % MSE = trace(C(:,:,1,1,1))/trace(Rscaled(:,:,1,1,1));
        
        %% Compute effective channel, effective channel estimation, and effective noise, depending on the type of network:
        switch net_type
            case 'CELLULAR'
                %% Generate combiner and precoding matrices for all UEs in each cell
                [~, W] =  functionCombinerPrecoder(Hhat_AP,B_AP,R_AP,batchsize,M,K,L,p, COMBINER, net_type);
                
                v_k = W(:,:,k,:);  %combining vectors for UE k and APs (normalize to make optimal s smaller) Size: M x batchsize x 1 x L
                hhat_k = Hhat_AP(:,:,k); %estimated channel between  UE k and all the APs. Size: LM x batchsize
                
                for batch_sample = 1:batchsize
                    %UL
                    H_AP_samp = squeeze(H_AP(:,batch_sample,:)); %pick out current sample, size: LM x K
                    v_k_samp = reshape(squeeze(v_k(:,batch_sample,:,:)),[L*M,1]); %pick out current sample and reshape. Size: LM x 1
                    
                    combinedChannels = (APassignment(:,k) .* v_k_samp)' * H_AP_samp;
                    g_ul((batch_indx-1)*batchsize + batch_sample) = combinedChannels(k); % Effective channel for UE k.
                    combinedChannels(k) = []; %remove the UEs channel to only have noise terms left
                    
                    ghat_ul((batch_indx-1)*batchsize +batch_sample) = (APassignment(:,k) .* v_k_samp)' * hhat_k(:,batch_sample); % Channel estimate after combining for UE k.
                    sigma_sq_ul((batch_indx-1)*batchsize + batch_sample) = vecnorm((APassignment(:,k) .* v_k_samp))^2 +  p*sum(abs(combinedChannels).^2); % Effective noise after combining.
                    
                    %DL
                    if L>1
                        W = reshape(permute(W,[1,4,2,3]),[L*M,batchsize,K]); % Size: LM x nbrOfRealizations x K
                    end
                    combinedChannels = H_AP(:,batch_sample,k)' * (APassignment.*squeeze(W(:,batch_sample,:))); % Precoder from cell managing each UE.
                    g_dl((batch_indx-1)*batchsize + batch_sample) = combinedChannels(k);
                    if UE_CSI == 2
                        ghat_dl((batch_indx-1)*batchsize + batch_sample) = functionChannelEstimatesDL(combinedChannels,length(batch_sample),np,pilotIndex,p,ESTIMATOR,k);
                    end
                    combinedChannels(k) = []; %remove the UEs channel to only have noise terms left
                    
                    sigma_sq_dl((batch_indx-1)*batchsize + batch_sample) = 1 +  p*sum(abs(combinedChannels).^2);
                end
                if UE_CSI ~= 2
                    ghat_dl((batch_indx-1)*batchsize + 1 : batch_indx*batchsize) = sum(conj(repmat(APassignment(:,k),[1,batchsize]) .*H_AP(:,:,k)) .* W(:,:,k), 1); % this will be used to estimate the channel hardening for DL
                end
            case 'DISTRIBUTED' 
                %% Generate combiner and precoding matrices for all UEs in each cell
                [~, W] =  functionCombinerPrecoder(Hhat_AP,B_AP,R_AP,batchsize,M,K,L,p, COMBINER, net_type);
                v_k = W(:,:,k,:);  %combining vectors for UE k and APs (normalize to make optimal s smaller) Size: M x batchsize x 1 x L
                for batch_sample = 1:batchsize
                    H_AP_samp = squeeze(H_AP(:,batch_sample,:)); %pick out current sample, size: LM x K
                    v_k_samp = reshape(v_k(:,batch_sample,:,:),[M,L]); %pick out current sample and reshape. Size: M x L
                    for l = 1:L
                        combinedChannels = v_k_samp(:,l)' * H_AP_samp(M*(l-1)+1:M*l,:); % Size L x K
                        g_ul_all_dist(l,(batch_indx-1)*batchsize + batch_sample,:) = combinedChannels; % Effective channels. Size: L x batch_sample x K
                    end
                    
                    % DL:
                    if L>1
                        W = reshape(permute(W,[1,4,2,3]),[L*M,batchsize,K]); % Size: LM x nbrOfRealizations x K
                    end
                    W_k_samp = squeeze(W(:,batch_sample,:)); %pick out current sample and reshape. Size: ML x K
                    combinedChannels = 1/sqrt(L) * H_AP(:,batch_sample,k)' * W_k_samp; % Size 1 x K
                    g_dl((batch_indx-1)*batchsize + batch_sample) = combinedChannels(k);
                    if UE_CSI == 2
                        ghat_dl((batch_indx-1)*batchsize + batch_sample) = functionChannelEstimatesDL(combinedChannels,length(batch_sample),np,pilotIndex,p,ESTIMATOR,k);
                    end
                    combinedChannels(k) = []; %remove the UEs channel to only have noise terms left
                    sigma_sq_dl((batch_indx-1)*batchsize + batch_sample) = 1 +  p*sum(abs(combinedChannels).^2);
                end
                if UE_CSI ~= 2
                    ghat_dl((batch_indx-1)*batchsize + 1 : batch_indx*batchsize) = sum(1/sqrt(L) * conj(H_AP(:,:,k)) .* W(:,:,k), 1); % this will be used to estimate the channel hardening for DL
                end
                % NOTE: The computation of g, ghat, and sigma is done after this for loop below since we need all the samples to compute hardening-based quantities 
            otherwise %'CENTRALIZED'
                %% Generate combiner and precoding matrices for all UEs in each cell
                [~, W] =  functionCombinerPrecoder(Hhat_AP,B_AP,R_AP,batchsize,M,K,L,p, COMBINER, net_type);
                
                v_k = W(:,:,k);  %combining vectors for UE k and APs (normalize to make optimal s smaller) Size: NL x batchsize x 1
                hhat_k = Hhat_AP(:,:,k); %estimated channel between  UE k and all the APs. Size: LN x batchsize x 1
                
                for batch_sample = 1:batchsize
                    %UL
                    H_AP_samp = squeeze(H_AP(:,batch_sample,:)); %pick out current sample, size: LN x K
                    v_k_samp = v_k(:,batch_sample); %pick out current sample and reshape. Size: LN x 1
                    
                    combinedChannels = v_k_samp' * H_AP_samp;
                    g_ul((batch_indx-1)*batchsize + batch_sample) = combinedChannels(k); % Effective channel for UE k.
                    
                    combinedChannels(k) = []; %remove the UEs channel to only have noise terms left
                    
                    ghat_ul((batch_indx-1)*batchsize +batch_sample) = v_k_samp' * hhat_k(:,batch_sample); % Channel estimate after combining for UE k.
                    sigma_sq_ul((batch_indx-1)*batchsize + batch_sample) = vecnorm(v_k_samp)^2 +  p*sum(abs(combinedChannels).^2); % Effective noise after combining.
                    
                    %DL
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
    end
 
    if strcmp(net_type,'DISTRIBUTED') % We need this part to implement the hardening combining at the CPU.
        % The function pagetimes will be used to perform matrix
        % multiplication for the K users at the same time. The same about
        % the function pagectranspose. 
        % The implementation of level 3 was started but never finished)
         
        %D_k = eye(L); % (for level 3) We use normalized combiners, so all the entries of this diagonal matrix are equal to 1.
        mean_g_ul_k = mean(g_ul_all_dist(:,:,k),2); %estimate channel hardening. Size: L x 1
        %mean_g_gh_ul = 1/nbrOfRealizations*pagemtimes(g_ul_all_dist,pagectranspose(g_ul_all_dist)); % (for level 3) Size: LxLxK ( Expectation inside inv in Eq. (21) "Making Cell-Free Massive MIMO Competitive With MMSE Processing and Centralized Implementation")
        % Concatenation of for loops to avoid problems with very large matrices:
        for batch_indx =1:batches
            for batch_sample = 1:batchsize
                % For level 3 of cooperation: 
                %a_k = (p*sum(mean_g_gh_ul,3) + D_k)\mean_g_ul_k; % Eq (21) in "Making Cell-Free Massive MIMO Competitive With MMSE Processing and Centralized Implementation"
                %a_k = a_k/(abs(sum(a_k))); % With this normalization works! But I'm not sure if it's correct!
                
                % For Level 2 of cooperation:
                a_k = 1/L*ones(L,1); % According to Eq. (23) and paragraph below in "Making Cell-Free Massive MIMO Competitive With MMSE Processing and Centralized Implementation"
                % The following is common in level 2 and 3:
                combinedChannels = squeeze(pagetimes(a_k',g_ul_all_dist(:,batch_sample,:))); % Size K x 1 (If >2021Ra, pagemtimes built-in function can be used instead)
                
                g_ul((batch_indx-1)*batchsize + batch_sample) = combinedChannels(k); % Effective channel for UE k.
                combinedChannels(k) = []; %remove the picked UE channel to only have noise terms left
                
                ghat_ul((batch_indx-1)*batchsize +batch_sample) = a_k'*mean_g_ul_k; % Channel estimate after combining for UE k.
                sigma_sq_ul((batch_indx-1)*batchsize + batch_sample) = vecnorm((a_k' * v_k_samp.').')^2 +  p*sum(abs(combinedChannels).^2); % Effective noise after combining.
                
            end

        end
    end
    
    if UE_CSI == 0
        ghat_dl = mean(ghat_dl) * ones(size(g_dl)); %estimate channel hardening DL
        % ghat_ul = mean(ghat_ul) * ones(size(g_ul)); %estimate channel hardening UL. To observe that the same is obtained in the UL if channel hardening is used as channel estimate after combining.
    end

time = toc;
disp(['Scenario + estimation + combiner takes: ' num2str(time) ' seconds']);

%-----------------------------------------------
% Estimate the uplink error probability
%-----------------------------------------------
tic

rate = infobits / n_ul;
f_ul = @(s) getErrorProbabilityUL(s, n_ul, p, rate, g_ul, ghat_ul, sigma_sq_ul);%only a tenth of samples used to find candidate
[eps_ul(pos_idx),  s_val_ul(pos_idx)] = golden_search(f_ul,1e-8,1,1e-8,eps_target);
time2 = toc;

disp(['UL search took: ' num2str(time2) ' seconds'])

%-----------------------------------------------
% Estimate the downlink error probability
%-----------------------------------------------
%if ~strcmp(net_type,'DISTRIBUTED') % The downlink for the cell-free distributed case is missing
rate = infobits / n_dl;

tic
f_dl = @(s) getErrorProbabilityDL(s,n_dl, p, rate, g_dl, ghat_dl, sigma_sq_dl); %only a tenth of samples used to find candidate
[eps_dl(pos_idx),  s_val_dl(pos_idx)] = golden_search(f_dl,1e-8,1,1e-8,eps_target);
time4=toc;
disp(['DL search took: ' num2str(time4) ' seconds'])
%end
disp(['UL: epsilon = ' num2str(eps_ul(pos_idx)) ' and s = ' num2str(s_val_ul(pos_idx))]);
disp(['DL: epsilon = ' num2str(eps_dl(pos_idx)) ' and s = ' num2str(s_val_dl(pos_idx))]);
fprintf('\n')

end

if DEBUG == 1 % We plot the figure. Add here any sanity check to be showed only in debug mode.
    
    keyboard
    figure
    [cdf, x]=ecdf(eps_ul);
    loglog(x,cdf); hold on;
    [cdf, x]=ecdf(eps_dl);
    loglog(x,cdf);
    
end

data.ASDdeg = ASDdeg;
data.L=L;
data.K=K;
data.M=M;
data.np=np;
data.n = n;
data.rate =  infobits/n_ul;
data.snr_db = rho_dBm;
data.nbrOfRealizations = nbrOfRealizations;
data.nbrOfPositions = nbrOfPositions;
data.COMBINER = COMBINER;
data.ESTIMATOR = ESTIMATOR;
data.net_type = net_type;
data.s_val_ul = s_val_ul;
data.avg_error_ul = eps_ul;
data.s_val_dl = s_val_dl;
data.avg_error_dl = eps_dl;
data.UE_CSI = UE_CSI;

%-------------------------------------------------
%       SAVE FILE
%-------------------------------------------------

filename = ['CDF_points_simNr_' num2str(simNr) '_' net_type '_' ESTIMATOR '_' COMBINER '_' antennaType '_M_' num2str(M) '_L_' num2str(L) '_K_' num2str(K) '_rho_' num2str(rho_dBm) '_np_' num2str(np) '_n_' num2str(n) '_UE_CSI_' num2str(UE_CSI) '_ASD_' num2str(ASDdeg) '.mat'];
save(filename, 'data', '-v7.3')

disp('Finished')

end

function avg_error = getErrorProbabilityUL(s, n, rho, rate, g_list, ghat_list, sigma_sq_list)
nbrOfRealizations = length(g_list);
eps_ul = nan(1, nbrOfRealizations);

b = rate*n;

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
    Is = preterm_ul + betaA_ul - betaB_ul;
    Vs = (betaA_ul - betaB_ul).^2 + 2*betaA_ul*betaB_ul*(1-nu_ul);
    
    eps_ul(i) = qfunc((n*Is-log(2^b-1)) ./ sqrt(n*Vs));
end
avg_error=mean(eps_ul);

end

function avg_error = getErrorProbabilityDL(s, n, rho, rate, g_list, ghat_list, sigma_sq_list)

b = rate*n;
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
    
    
    Is = preterm_dl + betaA_dl - betaB_dl;
    Vs = (betaA_dl - betaB_dl).^2 + 2*betaA_dl*betaB_dl*(1-nu_dl);
    
    eps_dl(i) = qfunc((n*Is-log(2^b-1)) ./ sqrt(n*Vs));
    
end

avg_error=  mean(eps_dl);

end
