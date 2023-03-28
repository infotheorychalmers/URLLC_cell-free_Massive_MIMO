function [R_AP,pilotIndex,gainOverNoisedB_AP,UEpositions,APpositions,pilotIndexDL] = generateSetup_hardening_study(L,K,N,np,ASDdeg,UE_POS,squareLength)
%Generate realizations of the simulation setup described in Section IV.
%
%This function was modified from the one developed as a part of the paper:
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
%L                  = Number of APs for the Cell-free system
%K                  = Number of UEs in the network
%N                  = Number of antennas per AP
%nbrOfSetups        = Number of setups with random UE locations
%
%OUTPUT:
%R_AP               = Matrix with dimension N x N x L x K x nbrOfSetups
%                     where (:,:,l,k,n) is the spatial correlation matrix
%                     between AP l and UE k in setup n, normalized by noise
%                     power
%pilotIndex         = Matrix with dimension K x nbrOfSetups containing the
%                     pilot assigned to the UEs in different setups
%gainOverNoisedB_AP = Matrix with dimension L x K x nbrOfSetups where
%                     element (l,k,n) is the channel gain (normalized by
%                     the noise power) between AP l and UE k in setup n
%UEpositions        = Vector of length K with UE positions, where the real
%                     part is the horizontal position and the imaginary
%                     part is the vertical position
%APpositions        = Vector of length L with the AP locations, measured in
%                     the same way as UEpositions

DEBUG = 0; % If DEBUG = 1, we plot the randomly generated setup. 

%% Define simulation setup

if nargin < 7
    %Size of the coverage area (as a square with wrap-around)
    squareLength = 150; %meter
elseif nargin < 6
    squareLength = 150; %meter
    UE_POS =7;
end
%Communication bandwidth
B = 20e6;

%Noise figure (in dB)
noiseFigure = 5; 

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Pathloss parameters
alpha = 36.7;
constantTerm = -30.5; % In original code

%Standard deviation of the shadow fading
sigma_sf = 4;

%Decorrelation distance of the shadow fading
decorr = 9;

%Height difference between an AP and a UE (In our Massive MIMO paper there
%is not vertical distance. However, here we don't have minimum distance. 
distanceVertical = 10;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = 1/2; %Half wavelength distance

%Number of APs per dimension on the grid
nbrAPsPerDim = sqrt(L);

%Distance between APs in vertical/horizontal direction
interAPDistance = squareLength/nbrAPsPerDim;

%Deploy APs on the grid
locationsGridHorizontal = repmat(interAPDistance/2:interAPDistance:squareLength-interAPDistance/2,[nbrAPsPerDim 1]);
locationsGridVertical = locationsGridHorizontal';
APpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);


%Prepare to save results
gainOverNoisedB_AP = zeros(L,K);
R_AP = zeros(N,N,L,K);
pilotIndex = zeros(1,K);

%Assign pilot sequences to UEs
pilotPattern = nan(1,K);
pilotSequences = 1:np;
for j = 1:K
    if isempty(pilotSequences)
        pilotSequences = 1:np; %create new available patterns
    end
    seqNr = datasample(pilotSequences,1,'Replace',false); %sample a pilot sequence
    pilotPattern(j) = seqNr; %assign UE j in cell i the sequence
    pilotSequences = pilotSequences(pilotSequences ~= seqNr); %remove pilot sequence from available sequences
end
pilotIndex = pilotPattern;

%Assign pilot sequences to APs (for DL. So far, orthogonal pilots for each AP)
pilotIndexDL = 1:L;

%Prepare to compute UE locations
%UEpositions = zeros(K,1);

%Prepare to store shadowing correlation matrix
shadowCorrMatrix = sigma_sf^2*ones(K,K);

% Different option for UE positions:
switch UE_POS
    % Good locations for cellular setting:
    case 1
        UEpositions = (0.75 + 1i*0.5) * squareLength; % User surrounded by APs with symmetry
    case 2
        UEpositions = (0.775 + 1i*0.5) * squareLength; % To induce a bit of assymetry wrt APs
    % Intermediate location for cellular setting:
    case 3
        UEpositions = (0.75 + 1i*0.75) * squareLength; % User surrounded by APs with symmetry
    case 4
        UEpositions = (0.775 + 1i*0.775) * squareLength; % To induce a bit of assymetry wrt APs
    % Bad locations for cellular setting
    case 5
        UEpositions = (0.5 + 1i*0.75) * squareLength; % User surrounded by APs with symmetry
    case 6
        UEpositions = (0.5 + 1i*0.775) * squareLength; % To induce a bit of assymetry wrt APs
    case 7
        UEpositions = (0.5 + 1i*0.5) * squareLength; % UE center of scenario. Always symmetery when deploying increasing # APs
    otherwise
        UEpositions = (rand + 1i*rand) * squareLength;
end

%Compute alternative AP locations by using wrap around
% wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
% wrapVertical = wrapHorizontal';
% wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
% APpositionsWrapped = repmat(APpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);

%Create correlated shadow fading realizations
shadowAPrealizations = sqrtm(shadowCorrMatrix)*randn(K,L);

%Go through all UEs
for k = 1:K

    %Compute distances assuming for the UE to the APs
%     [distanceAPstoUE,whichpos] = min(abs(APpositionsWrapped -
%     repmat(UEpositions(k),size(APpositionsWrapped))),[],2); % with wrap-around topology
    distanceAPstoUE = abs(APpositions - UEpositions(k));

    distances = sqrt(distanceVertical^2+distanceAPstoUE.^2);

    %Compute the channel gain divided by the noise power (in dB)
    %gainOverNoisedB_AP(:,k,n) = constantTerm - alpha*log10(distances) + shadowAPrealizations(k,:)' - noiseVariancedBm; % with shadowing
    gainOverNoisedB_AP(:,k) = constantTerm - alpha*log10(distances) - noiseVariancedBm; % no shadowing.

    %Go through all APs
    for l = 1:L

        %Compute nominal angle between UE k and AP l
%         angletoUE = angle(UEpositions(k)-APpositionsWrapped(l,whichpos(l))); % with wrap-around topology
        angletoUE = angle(UEpositions(k)-APpositions(l));

        %Generate normalized spatial correlation matrix using the local scattering model
        R_AP(:,:,l,k) = db2pow(gainOverNoisedB_AP(l,k))*functionRlocalscattering(N,angletoUE,ASDdeg,antennaSpacing);

    end

end

if DEBUG
    figure
    scatter(real(APpositions(:)),imag(APpositions(:)),10,'b','filled');
    xticks(0:squareLength/nbrAPsPerDim:squareLength)
    yticks(0:squareLength/nbrAPsPerDim:squareLength)
    grid on;
    hold on
    plot(real(UEpositions), imag(UEpositions),'*r')
end

