function [ghat_dl] = functionChannelEstimatesDL(combinedChannels,nbrOfRealizations,np,pilotIndex,p,ESTIMATOR,UE_idx)
%Generate the channel realizations and estimates of these channels for all
%UEs in the entire network. The channels are modeled as correlated
%Rayleigh fading and the MMSE estimator is used.
%
%This function was developed as a part of the paper:
%
%Emil Bjornson, Luca Sanguinetti, "Making Cell-Free Massive MIMO
%Competitive With MMSE Processing and Centralized Implementation,"
%IEEE Transactions on Wireless Communications, To appear.
%
%Download article: https://arxiv.org/abs/1903.10611
%
%This is version 1.0 (Last edited: 2019-03-19)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%paper as described above.
%
%INPUT:
%combinedChannels  = Effective DL channels after precoding 1 x K
%nbrOfRealizations = Number of channel realizations
%M                 = Number of antennas per AP
%np                = Number of orthogonal pilots
%pilotIndex        = Vector containing the pilot assigned to each UE
%p                 = Uplink transmit power per UE (same for everyone)
%ESTIMATOR         = Channel estimation technique [MMSE,LS]
% 
%OUTPUT:
%ghat_dl      = Matrix with dimension nbrOfRealizations x 1            


%% Perform channel estimation
K = length(pilotIndex);
%Generate realizations of normalized noise
Np = sqrt(0.5)*(randn(nbrOfRealizations,K) + 1i*randn(nbrOfRealizations,K));

%Prepare to store results
ghat_dl = zeros(nbrOfRealizations,1);
   
%Compute processed pilot signal for analyzed UE that use pilot t
t = pilotIndex(UE_idx);
yp = sqrt(p)*np*sum(combinedChannels(find(t==pilotIndex))) + sqrt(np)*Np(:,UE_idx);

%Compute the matrix that is inverted in the MMSE estimator
PsiInv = (p*np*sum(t==pilotIndex) + 1);

if strcmp(ESTIMATOR, 'MMSE')%Compute the MMSE estimate
%     RPsi = 1 / PsiInv;
%     ghat_dl = sqrt(p)*RPsi*yp; This one is wrong, so we always use LS
%     with DL pilots. 
    A_LS = 1/(sqrt(p)*np);
    ghat_dl = A_LS*yp;
elseif strcmp(ESTIMATOR, 'LS')
    A_LS = 1/(sqrt(p)*np);
    ghat_dl = A_LS*yp;
end

 
end