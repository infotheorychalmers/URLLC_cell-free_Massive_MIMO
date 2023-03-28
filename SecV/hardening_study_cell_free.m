function data = hardening_study_cell_free(L_tot, K, M, infobits, n, np_ul,np_dl, rho_dBm, ESTIMATOR, COMBINER, nbrOfRealizations, antennaType, net_type, simNr_vec, UE_CSI, CSI_AP,UE_POS)
vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix

DEBUG = 0;
rng('shuffle')

if DEBUG == 1
    %FLAGS

    K = 1; % Number of UEs in the network
    net_type = 'CENTRALIZED'; % type of network [CENTRALIZED,DISTRIBUTED]. Only centralized here.
    L_tot = 144;
    M = 4;
    np_ul = 40; %number of UL pilots
    np_dl = 40; %number of DL pilots
    n = 300; %total number of channle uses
    infobits = 8*20; %transmission rate
    rho_dBm = -10; %transmit power [dBm]
    nbrOfRealizations=1e5; %number of saddlepoint realizations
    UE_POS = 7; % position of the UE (scaled by squareLength) - 1:(0.75 + 1i*0.5); 2:(0.775 + 1i*0.5); 3:(0.75 + 1i*0.75); 4:(0.775 + 1i*0.775); 5:(0.5 + 1i*0.75); 6:(0.5 + 1i*0.775)
    COMBINER = 'MR'; %what combiner to use [MR, M-MMSE, RZF]
    ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
    antennaType = 'ULA';

    simNr_vec = 1:10;
    UE_CSI=2; % UE_CSI=0 means decoding relying on channel hardening
    % UE_CSI=1 means same CSI at the BS and UE
    % UE_CSI=2 means pilots in the DL
    CSI_AP = 0;
end
L_list = (2:2:sqrt(L_tot)).^2;

eps_target = 1e-5;
squareLength = 150; % size of scenario (meter)
ASDdeg = 25; %Angular standard deviation around the nominal angle (measured in degrees)

p = 10^(rho_dBm/10);

%% Run Simulation

% Initializations:
g_ul = nan(nbrOfRealizations,1);
g_dl = g_ul;
ghat_dl = g_ul;
sigma_sq_dl = g_ul;

%% Generate setup and compute channel statistics:

% Further initializations:
hardening = nan(1,length(L_list));
eps_dl = nan(1,length(L_list));
s_val_dl = nan(1,length(L_list));
APpositions_scalable = nan(L_tot,length(L_list));
for ss=1:length(simNr_vec)
    simNr = simNr_vec(ss);
    for ll = 1:length(L_list)
        tic
        L = L_list(ll);
        if UE_CSI ~= 2
            np_dl = 0;
        end
        n_ul = n/2 - np_ul; %length of UL phase
        n_dl = n/2 - np_dl; %length of DL phase

        [R_AP,pilotIndex,~,UEpositions,APpositions,pilotIndexDL] = generateSetup_hardening_study(L,K,M,np_ul,ASDdeg,UE_POS,squareLength);

        if DEBUG
            %         figure
            %         scatter(real(APpositions),imag(APpositions),10,'b','filled');
            %         %             xticks(0:squareLength/nbrAPsPerDim:squareLength)
            %         %             yticks(0:squareLength/nbrAPsPerDim:squareLength)
            %         axis([0 squareLength 0 squareLength])
            %         grid on;
            %         hold on
            %         plot(real(UEpositions), imag(UEpositions),'*r')
        end
        %R_values = squeeze(R_AP) % Uncomment to display the correlation matrices
        k = randi(K,1); %pick a random UE
        % Create batches:
        batchsize = min(1e2, nbrOfRealizations);
        batches = nbrOfRealizations / batchsize;

        %Need to do this for memory issues
        for batch_indx =1:batches

            %% Generate channel and channel estimates
            [Hhat_AP,H_AP,B_AP] = functionChannelEstimates(R_AP,batchsize,L,K,M,np_ul,pilotIndex,p,ESTIMATOR);
            % MSE = trace(C(:,:,1,1,1))/trace(Rscaled(:,:,1,1,1));
            %hardening_prev_combining (ll,mm) = var(vecnorm(H_AP).^2)/mean(vecnorm(H_AP).^2)^2;
            switch net_type
                case 'CENTRALIZED'
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
                            ghat_dl((batch_indx-1)*batchsize + batch_sample) = functionChannelEstimatesDL(combinedChannels,length(batch_sample),np_dl,pilotIndex,p,ESTIMATOR,k);
                        end
                        combinedChannels(k) = []; %remove the UEs channel to only have noise terms left

                        sigma_sq_dl((batch_indx-1)*batchsize + batch_sample) = 1 +  p*sum(abs(combinedChannels).^2);
                    end
                    if UE_CSI ~= 2
                        ghat_dl((batch_indx-1)*batchsize + 1 : batch_indx*batchsize) = sum(conj(H_AP(:,:,k)) .* W(:,:,k), 1); % this will be used to estimate the channel hardening for DL
                    end
                otherwise  %'DISTRIBUTED'
                    if CSI_AP == 0
                        [~, W] =  functionCombinerPrecoder(Hhat_AP,B_AP,R_AP,batchsize,M,K,L,p, COMBINER, net_type);
                    else
                        [~, W] =  functionCombinerPrecoder(H_AP,B_AP,R_AP,batchsize,M,K,L,p, COMBINER, net_type);
                    end
                    for batch_sample = 1:batchsize
                        W_k_samp = reshape(permute(squeeze(W(:,batch_sample,:,:)),[1,3,2]),[M*L,K]); %pick out current sample and reshape. Size: ML x K
                        combinedChannels = 1/sqrt(L)*H_AP(:,batch_sample,k)' * W_k_samp; % Size 1 x K
                        %                 g_ul_all_distDL((batch_indx-1)*batchsize + batch_sample,:) = combinedChannels; % Effective channels. Size: batch_sample x K

                        g_dl((batch_indx-1)*batchsize + batch_sample) = combinedChannels(k);
                        if UE_CSI == 2
                            ghat_dl((batch_indx-1)*batchsize + batch_sample) = functionChannelEstimatesDL(combinedChannels,length(batch_sample),np_dl,pilotIndex,p,ESTIMATOR,k);
                        end
                        combinedChannels(k) = []; %remove the UEs channel to only have noise terms left

                        sigma_sq_dl((batch_indx-1)*batchsize + batch_sample) = 1 +  p*sum(abs(combinedChannels).^2);
                    end
                    if UE_CSI ~= 2
                        if L>1
                            W = reshape(permute(W,[1,4,2,3]),[L*M,batchsize,K]); % Size: LM x nbrOfRealizations x K
                        end
                        ghat_dl((batch_indx-1)*batchsize + 1 : batch_indx*batchsize) = sum(1/sqrt(L)*conj(H_AP(:,:,k)) .* W(:,:,k), 1); % this will be used to estimate the channel hardening for DL
                    end

            end



        end



        %check channel hardening
        hardening(ll) = var(abs(ghat_dl).^2)/mean(abs(ghat_dl).^2)^2;
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
        [eps_dl(ll),  s_val_dl(ll)] = golden_search(f_dl,1e-8,1,1e-8,eps_target);
        time4=toc;
        disp(['DL search took: ' num2str(time4) ' seconds'])

        %disp(['UL: epsilon = ' num2str(eps_ul(ll)) ' and s = ' num2str(s_val_ul(ll))]);
        disp(['DL: epsilon = ' num2str(eps_dl(ll)) ' and s = ' num2str(s_val_dl(ll))]);
        fprintf('\n')
        if eps_dl(ll) < eps_target/10
            break
        end
    end


    %-------------------------------------------------
    %       FIGURES
    %-------------------------------------------------

    % Plot variance channel hardening after combining
    figure
    plot(M*L_list, hardening)
    title(['centralized cell-free network with M = ' num2str(M) ', AP CSI = ' num2str(CSI_AP) ', UE CSI = ' num2str(UE_CSI)])
    xlabel('total number of antennas in the system, LM')
    ylabel('Normalized variance hardening metric')
    grid on
    % Plot downlink average error probability
    figure
    semilogy(M*L_list,eps_dl)
    ylim([1e-5,1])
    title(['centralized cell-free network with M = ' num2str(M) ', AP CSI = ' num2str(CSI_AP) ', UE CSI = ' num2str(UE_CSI)])
    xlabel('total number of antennas in the system, LM')
    ylabel('downlink average error probability')
    grid on


    time = toc;
    disp(['Scenario + estimation + combiner takes: ' num2str(time) ' seconds']);


    %-------------------------------------------------
    %       SAVE FILE
    %-------------------------------------------------

    data.ASDdeg = ASDdeg;
    data.L_tot = L_tot;
    data.L_list = L_list(~isnan(eps_dl));
    data.K=K;
    data.M=M;
    data.np_ul=np_ul;
    data.np_dl=np_dl;
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
    data.APpositions = APpositions_scalable(:,1:length(data.L_list));
    data.squareLength = squareLength;
    % %-------------------------------------------------
    % %       SAVE FILE
    % %-------------------------------------------------

    filename = ['hardening_simNr_' num2str(simNr) '_cell-free_' net_type '_' ESTIMATOR '_' COMBINER '_' antennaType '_M_' num2str(M) '_L_' num2str(L_tot) '_CSI_AP_' num2str(CSI_AP) '_CSI_UE_' num2str(UE_CSI) '_K_' num2str(K) '_UE_POS_' num2str(UE_POS) '_rho_' num2str(rho_dBm) '_np_ul_' num2str(np_ul) '_np_dl_' num2str(np_dl) '_n_' num2str(n) '.mat'];

    save(filename, 'data', '-v7.3')

    disp('Finished')

end
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

