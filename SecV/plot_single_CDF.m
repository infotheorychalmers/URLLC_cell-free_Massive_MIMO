clear all
%FLAGS
np = [40]; %number of pilots
n = 300; %total number of channle uses
b = 8*20; %transmission rate
snr_db = 10; %transmit power [dBm] (10,15,20)
nbrBSs = 1; % Number of BSs in the cellular network. Cannot be changed!
%Msyst = 400;
%M = Msyst/nbrBSs; % total number of antennas managed by each BS
L = 4;
N = 25;
M = N*L;
K=40; % number of users in the entire scenario
COMBINER = 'M-MMSE'; %what combiner to use [MR, M-MMSE,RZF]
ESTIMATOR = 'MMSE'; %what combiner to use [LS, MMSE]
antennaType = 'ULA';
ASDdeg=25; % (10,20,30)
level = 1;


filename = [num2str(level) '_' ESTIMATOR '_' COMBINER '_' antennaType '_M_' num2str(M) '_N_' num2str(N) '_K_' num2str(K) '_L_' num2str(L) '_SNR_' num2str(snr_db) '_np_' num2str(np) '_n_' num2str(n) '_ASD_' num2str(ASDdeg) '.mat'];
%filename2 = 'MMSE_RZF_M_25_SNR_10_np_40_n_300_ASD_25_K_10_L_4';
folder = [pwd '/data/'];

simNr = 1:4;

for i = 1:length(simNr)
    cur_file = ['simNr_' num2str(simNr(i)) '_Multicell_cellfree_CDF_level_' filename];
    %cur_file = ['simNr_' num2str(simNr(i)) '_Multicell_ULA_CDF_' filename2];
    
    try
        load([folder, cur_file])
    catch
        disp(['The file ' cur_file ' was not found so we skip it.'])
        continue
    end
    
    Nl = length(data.avg_error_ul);
    
    eps_ul((i-1)*Nl+1 : i*Nl) = data.avg_error_ul;
    %             eps_dl((i-1)*Nl+1 : i*Nl) = data.avg_error_dl;
end

%Obtain UL params
[cdf_ul, x_ul] = ecdf(eps_ul);

figure
loglog(x_ul, cdf_ul); hold on;

savename_ul = ['CDF_Multicell_UL_cellfree_level_' filename '.csv']
csvwrite(savename_ul, [x_ul, cdf_ul]);




