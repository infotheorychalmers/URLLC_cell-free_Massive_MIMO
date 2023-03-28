clear all
%FLAGS
K = 40; % Number of UEs in the network
n = 300; %total number of channel uses
b = 8*20; %transmission rate
rho_dBm_list = [-20:2: 10 20 30 34]; %transmit power [dB]
np = K;
L = 25; % Number of access points in the network
M = 4; %Total number of antennas per AP 
COMBINER = 'M-MMSE'; %what combiner to use [MR, M-MMSE,RZF]
ESTIMATOR = 'LS'; %what combiner to use [LS, MMSE]
antennaType = 'ULA';
net_type = 'CENTRALIZED';

simNr = 1;
UE_CSI =2;

ASD = 25;
eps_target = 1e-5;


ul_availability = nan(length(rho_dBm_list), 1);
dl_availability = nan(length(rho_dBm_list), 1);


for rhoindx = 1:length(rho_dBm_list)
    rho_dBm=rho_dBm_list(rhoindx);
    
    folder = [pwd '/data_rho/'];

    for i = 1:length(simNr)
        filename = ['CDF_points_simNr_' num2str(simNr(i)) '_' net_type '_' ESTIMATOR '_' COMBINER '_' antennaType '_M_' num2str(M) '_L_' num2str(L) '_K_' num2str(K) '_rho_' num2str(rho_dBm) '_np_' num2str(np) '_n_' num2str(n) '_UE_CSI_' num2str(UE_CSI) '_ASD_' num2str(ASD)];
        cur_file = [filename '.mat'];
        try
            load([folder, cur_file])
        catch
            disp(['The file ' cur_file ' was not found so we skip it.'])
            continue
        end
        
        Nl = length(data.avg_error_ul);
        
        eps_ul((i-1)*Nl+1 : i*Nl) = data.avg_error_ul;
        eps_dl((i-1)*Nl+1 : i*Nl) = data.avg_error_dl;
    end
    
    %Obtain UL params
    [cdf_ul, x_ul] = ecdf(eps_ul);
    [cdf_dl, x_dl] = ecdf(eps_dl);
    
    %                     figure(Mindx)
    %                     loglog(x_ul, cdf_ul); hold on;
    %                     loglog(x_dl, cdf_dl);
    
    [~ , target_indx] = min(abs(x_ul-eps_target));
    
    if max(x_ul) < eps_target
        target_indx = length(cdf_ul);
    end
    
    ul_availability(rhoindx) = cdf_ul(target_indx);
    
    [~ , target_indx] = min(abs(x_dl-eps_target));
    
    if max(x_dl) < eps_target
        target_indx = length(cdf_dl);
    end
    
    dl_availability(rhoindx) = cdf_dl(target_indx);
    
end

    figure;hold on
    plot(rho_dBm_list, ul_availability);hold on;
    plot(rho_dBm_list, dl_availability);
    grid on;
    savename_ul = ['UL_' filename '.csv']
    savename_dl = ['DL_' filename '.csv']
%     csvwrite(savename_ul, [rho_dBm_list' , ul_availability]);
%     csvwrite(savename_dl, [rho_dBm_list' , dl_availability]);



