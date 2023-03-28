clear all
%FLAGS
K = 1; % Number of UEs in the network
net_type = 'CELLULAR'; % type of network [CELLULAR,CENTRALIZED]. Only cellular here.
L = 1; % Number of APs in the small-cell network
M_tot = 144; % Number of antennas per AP
np = 40; %number of pilots
n = 300; %total number of channle uses
infobits = 8*20; %transmission rate
rho_dBm = -10; %transmit power [dBm]
UE_POS_list = 1:6; % position of the UE (scaled by squareLength) - 1:(0.75 + 1i*0.5); 2:(0.775 + 1i*0.5); 3:(0.75 + 1i*0.75); 4:(0.775 + 1i*0.775); 5:(0.5 + 1i*0.75); 6:(0.5 + 1i*0.775)
COMBINER = 'MR'; %what combiner to use [MR, M-MMSE, RZF]
ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
antennaType = 'ULA';

M_list = [2:sqrt(M_tot)].^2;

simNr = 1;
UE_CSI=2; % UE_CSI=1 means same CSI at the BS and UE
CSI_AP = 0;

eps_target = 1e-5;
squareLength = 150; % size of scenario (meter)
ASDdeg = 25; %Angular standard deviation around the nominal angle (measured in degrees)
n_ul = (n - np) / 2; %length of UL phase
n_dl = n - np - n_ul; %length of DL phase
p = 10^(rho_dBm/10);

%Number of APs per dimension on the grid
nbrAPsPerDim = 1;

%Distance between APs in vertical/horizontal direction
interAPDistance = squareLength/nbrAPsPerDim;

%Deploy APs on the grid
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
% Save positions of all APs:
savename_BS = ['BSposition_hardening_cellular.csv']
csvwrite(savename_BS, [real(APpositions),imag(APpositions)]);
for UE_POS = 1:length(UE_POS_list)   
    filename = ['hardening_simNr_' num2str(simNr) '_' net_type '_' ESTIMATOR '_' COMBINER '_' antennaType '_M_' num2str(M_tot) '_L_' num2str(L) '_CSI_AP_' num2str(CSI_AP) '_CSI_UE_' num2str(UE_CSI) '_K_' num2str(K) '_UE_POS_' num2str(UE_POS) '_rho_' num2str(rho_dBm) '_np_' num2str(np) '_n_' num2str(n)];
    cur_file = [filename '.mat'];
    folder = [pwd '/data/'];
    
    try
        load([folder, cur_file])
    catch
        disp(['The file ' cur_file ' was not found so we skip it.'])
        continue
    end
    
    % Plot BS
    figure
    scatter(real(APpositions),imag(APpositions),10,'b','filled');
    hold on
    axis([0 squareLength 0 squareLength])
    grid on;
    % Save UE position
    savename_UE = ['UEposition_' filename '.csv'];
    csvwrite(savename_UE, [real(data.UEpositions),imag(data.UEpositions)]);
    % Plot position of UE
    plot(real(data.UEpositions),imag(data.UEpositions),'*r','MarkerSize',20)
    
    % Plot hardening metric and DL average error as a function of # active APs
    figure
%     plot(data.M_list, data.hardening);
    plot(data.M_list, data.hardening);
    figure
    semilogy(data.M_list, data.avg_error_dl);
    % Save hardening metric and DL average error as a function of # active APs
    savename_hard = [filename '.csv'];
    csvwrite(savename_hard, [data.M_list',data.hardening']);
    savename_err = ['DL_avg_err_' filename '.csv'];
    csvwrite(savename_err, [data.M_list',data.avg_error_dl']);
end




