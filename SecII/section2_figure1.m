%Empty workspace and close figures
close all;
clear;

%Set the side length of the simulation area
squareLength = 150;

%Total number of antennas in all setups
nbrOfAntennas = 64;

%Number of UEs in the simulation setup
K_values = [1,2:2:16 20:5:40];

%Number of realizations of the random UE locations
nbrOfSetups = 5000;

%Combiner
combiner = 'MR';

%%  Transmit parameters
%transmit power [dBm]
pow_dBm = -10;

%noise power [dBm]
sigma2_dBm = -96;

%Codeword length
codewordLength = 100;

%number of bits
nbrOfbits = 3*20;

%number of samples in UL
nbrOfSamplesUl = 100;

%number of bits
rate = nbrOfbits/codewordLength;

%Target error probability
eps_target = 1e-5;

%%  Network layout
%Set the AP location of the cellular Massive MIMO setup
APcellular = squareLength/2 + 1i*squareLength/2;

%Set the AP locations for the small-cell and cell-free setups
APperdim = sqrt(nbrOfAntennas);
APcellfree = linspace(squareLength/APperdim,squareLength,APperdim)-squareLength/APperdim/2;
APcellfree = repmat(APcellfree,[APperdim 1]) + 1i*repmat(APcellfree,[APperdim 1])';

% %%Plot network layout
% figure
% plot(real(APcellular),imag(APcellular),'rd','MarkerFaceColor','r',...
%     'MarkerSize',10);  hold on;
% plot(real(APcellfree),imag(APcellfree),'ko','MarkerFaceColor','k',...
%     'MarkerSize',8);
% legend({'Macro BS','APs'},'Interpreter','latex','Location','NorthEastOutside');
% set(gca,'fontsize',16);

%%  SNR computation
%Define a function to compute the SNR as function of the horizontal distance
%measured in meter. The AP is 10 meter above the UE.
SNR = @(hor_dist) db2pow(pow_dBm - sigma2_dBm - 30.5 - 36.7*log10(sqrt(hor_dist.^2+10^2)));

%% Prepare to store simulation results
eps_mMIMO    = zeros(nbrOfSetups,length(K_values));
eps_cellfree = zeros(nbrOfSetups,length(K_values));
eps_small    = zeros(nbrOfSetups,length(K_values));

availability_mMIMO = zeros(length(K_values),1);
availability_cellfree = zeros(length(K_values),1);
availability_small = zeros(length(K_values),1);

%% Go through all numbers of UEs
for ii = 1:length(K_values)
    
    %Output simulation progress
    disp([num2str(ii) ' UEs out of ' num2str(length(K_values))]);
    
    % Number of UEs
    K = K_values(ii);
    
    %Generate the random UE locations for all setups
    UElocations = squareLength*(rand(nbrOfSetups,K)+1i*rand(nbrOfSetups,K));
    
    %% Go through all random realizations of the UE locations
    for n = 1:nbrOfSetups
        
        if mod(n-1,100) == 0
            %Output simulation progress
            disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        end
        
        %Prepare to save channel matrix in the cellular Massive MIMO and
        %cell-free setups
        channelCellular = zeros(nbrOfAntennas,K);
        channelCellfree = zeros(nbrOfAntennas,K);
        
        for k = 1:K
            
            %Cellular Massive MIMO
            distanceCellular = abs(APcellular - UElocations(n,k));
%             channelCellular(:,k) = sqrt(SNR(distanceCellular));% Deterministic
            channelCellular(:,k) = sqrt(SNR(distanceCellular))*exp(1i*2*pi*rand(nbrOfAntennas,1));% Unitary            
%             channelCellular(:,k) = sqrt(0.5*SNR(distanceCellular)).*(randn(nbrOfAntennas,1)+1i*randn(nbrOfAntennas,1));% Gaussian
            %Cell-free Massive MIMO
            distanceCellfree = abs(APcellfree(:) - UElocations(n,k));
%             channelCellfree(:,k) = sqrt(SNR(distanceCellfree));
            channelCellfree(:,k) = sqrt(SNR(distanceCellfree)).*exp(1i*2*pi*rand(nbrOfAntennas,1));% Unitary
%             channelCellfree(:,k) = sqrt(0.5*SNR(distanceCellfree)).*(randn(nbrOfAntennas,1)+1i*randn(nbrOfAntennas,1));% Gaussian

            
        end
        
        %Compute error probability with mMIMO
        [eps_mMIMO(n,ii),  s_val_mMIMO] = functionComputeErrProb(channelCellular, ...
            nbrOfSamplesUl, rate, eps_target,combiner);
        
        %Compute error probability with cell-free
        [eps_cellfree(n,ii), s_val_cellfree] = functionComputeErrProb(channelCellfree, ...
            nbrOfSamplesUl, rate, eps_target,combiner);
        
        eps = zeros(nbrOfAntennas,1);
        
        %Compute error probability with small-cell
        for m = 1:nbrOfAntennas
            
            [eps(m), s_val] = functionComputeErrProb(channelCellfree(m,:), ...
                nbrOfSamplesUl, rate, eps_target,combiner);
            
        end
        
        %Assign each UE to the AP providing the smallest error probability
        eps_small(n,ii) = min(eps,[],1);
        
    end
    
    % Compute network availability for mMIMO
    availability_mMIMO(ii) = functionComputeAvailability(eps_mMIMO(:,ii), eps_target);
    % Compute network availability for cell-free
    availability_cellfree(ii) = functionComputeAvailability(eps_cellfree(:,ii), eps_target);
    % Compute network availability for small-cells
    availability_small(ii) = functionComputeAvailability(eps_small(:,ii), eps_target);
    
end

%% Plot simulation results
fig1 = figure(1);
hold on; box on; grid on;
plot(K_values,availability_cellfree,'b--','LineWidth',2);
plot(K_values,availability_small,'r-.','LineWidth',2);
plot(K_values,availability_mMIMO,'k','LineWidth',2);
legend({'Cell-free','Small cells','Massive MIMO'},'Interpreter','latex','Location','SouthEast');
set(gca,'fontsize',18);
xlim([min(K_values), max(K_values)]);
ylim([0., 1]);
xlabel('Number of UEs, $K$','Interpreter','latex');
ylabel('Network Availability, $\eta$','Interpreter','latex');
fig1.Position(3:4) = [550 350];

% save

savename = ['UL_simple_channel_K_vs_net_avail_' combiner 'm10dBm.csv']
csvwrite(savename, [K_values' , availability_cellfree, availability_small, availability_mMIMO]);