clear all
%FLAGS
K = 1; % Number of UEs in the network
net_type = 'CENTRALIZED'; % type of network [DISTRIBUTED,CENTRALIZED]. Only centralized here.
L_tot = 144;
M = 1;
np_ul = 40; %number of UL pilots
np_dl = 0; %number of DL pilots
n = 300; %total number of channle uses
infobits = 8*20; %transmission rate
rho_dBm = -10; %transmit power [dBm]
nbrOfRealizations=1e5; %number of saddlepoint realizations
UE_POS_list = 1:7; % position of the UE (scaled by squareLength) - 1:(0.75 + 1i*0.5); 2:(0.775 + 1i*0.5); 3:(0.75 + 1i*0.75); 4:(0.775 + 1i*0.775); 5:(0.5 + 1i*0.75); 6:(0.5 + 1i*0.775)
COMBINER = 'MR'; %what combiner to use [MR, M-MMSE, RZF]
ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
antennaType = 'ULA';

simNr = 1;
UE_CSI=0; % UE_CSI=1 means same CSI at the BS and UE
CSI_AP = 0;

eps_target = 1e-5;
squareLength = 150; % size of scenario (meter)
ASDdeg = 25; %Angular standard deviation around the nominal angle (measured in degrees)
if UE_CSI ~= 0
    n_ul = (n - np_ul-np_dl) / 2; %length of UL phase
    n_dl = n - np_ul - np_dl - n_ul; %length of DL phase
else
    n_ul = (n - np_ul) / 2; %length of UL phase
    n_dl = n - np_ul - n_ul; %length of DL phase
end

p = 10^(rho_dBm/10);

for UE_POS = 1:length(UE_POS_list)
    
    filename = ['hardening_simNr_' num2str(simNr) '_cell-free_' net_type '_' ESTIMATOR '_' COMBINER '_' antennaType '_M_' num2str(M) '_L_' num2str(L_tot) '_CSI_AP_' num2str(CSI_AP) '_CSI_UE_' num2str(UE_CSI) '_K_' num2str(K) '_UE_POS_' num2str(UE_POS) '_rho_' num2str(rho_dBm) '_np_ul_' num2str(np_ul) '_np_dl_' num2str(np_dl) '_n_' num2str(n)];
    cur_file = [filename '.mat'];
    folder = [pwd '/data/'];
    
    try
        load([folder, cur_file])
    catch
        disp(['The file ' cur_file ' was not found so we skip it.'])
        continue
    end
    
    % Save UE position
    savename_UE = ['UEposition_' filename '.csv'];
    csvwrite(savename_UE, [real(data.UEpositions),imag(data.UEpositions)]);
    
    if UE_POS == 1 % The AP positions are the same for all UE position. We only save it once. 
        for ll = 1:size(data.L_list,2)
            
            L = data.L_list(ll);
            %Number of APs per dimension on the grid
            nbrAPsPerDim = sqrt(L);
            
            %Distance between APs in vertical/horizontal direction
            interAPDistance = squareLength/nbrAPsPerDim;
            
            %Deploy APs on the grid
            locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
            locationsGridVertical = locationsGridHorizontal';
            APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
            % Plot positions of all APs in the network
            figure;
            scatter(real(APpositions),imag(APpositions),10,'b','filled');
            hold on
            axis([0 squareLength 0 squareLength])
            grid on;
            % Plot position of UE
            plot(real(data.UEpositions),imag(data.UEpositions),'*r','MarkerSize',20)
            % Save positions of APs
            savename_APs = ['APpositions_cell-free_dens' num2str(ll) '_' filename '.csv'];
            csvwrite(savename_APs, [real(APpositions),imag(APpositions)]);
        end
    end
    if UE_CSI == 0
    % Plot hardening metric and DL average error as a function of # active APs
    figure
    plot(data.L_list, data.hardening);
    end
    figure
    semilogy(data.L_list, data.avg_error_dl);
    % Save hardening metric and DL average error as a function of # active APs
    if UE_CSI == 0
    savename_hard = [filename '.csv'];
    csvwrite(savename_hard, [data.M*data.L_list',data.hardening']);
    end
    savename_err = ['DL_avg_err_' filename '.csv'];
    csvwrite(savename_err, [data.M*data.L_list',data.avg_error_dl']);
end




