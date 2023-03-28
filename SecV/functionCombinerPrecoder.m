function [V,W] = functionCombinerPrecoder(Hhat,B,R,nbrOfRealizations,N,K,L,p, combiner, net_type)
%Compute UL SE for different receive combining schemes using Theorem 4.1.
%
%INPUT:
%Hhat              = Matrix with dimension N*L x nbrOfRealizations x K
%                    where (:,n,k) is the estimated collective channel from
%                    all BSs to UE k at channel realization n.
%R                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix between BS
%                    l and UE k in setup n, normalized by the noise power
%B                 = Matrix with dimension N x N x L x K where
%                    (:,:,l,k) is the spatial correlation matrix of the
%                    estimate between BS l and UE k in setup n, normalized
%                    by the noise power
%nbrOfRealizations = Number of channel realizations
%N                 = Number of antennas per AP
%K                 = Total number of UEs
%L                 = Number of APs
%p                 = Matrix K x 1 where element k is the uplink transmit
%                    power of UE k (If it is a scalar, the same value is
%                    used for all users)
%combiner          = ['MR', 'RZF', 'MMSE']
%net_type          = % type of network [CELLULAR,DISTRIBUTED,CENTRALIZED] 
%
%OUTPUT:
%V                 = if level=1, matrix of dimension N x nbrOfRealizations x K x L containing all
%                   the combining vectors for each realization of user k wrt to AP L.
%                   if level=4, matrix of dimension NL x nbrOfRealizations x K containing
%                   the combining vectors for each realization of user k for the fully centralized scenario.
%W                 = recoding matrix. Same as V but normalized by the vector norm.

vecnorm = @(A)  sqrt(sum(A.*conj(A),1)); %norm of each column in matrix


%Store identity matrices of different sizes
eyeN = eye(N);
eyeK = eye(K);
eyeNL = eye(N*L);

%Compute sum of all estimation error correlation matrices at every BS
C_tot = zeros(N,N,L);
%C1_tot = zeros(N,N,L);
for k = 1:K
    C_tot = C_tot + p(1)*(R(:,:,:,k)-B(:,:,:,k)); % p(k) if users use different powers
    %C1_tot = C1_tot + p1(k)*(R(:,:,:,k)-B(:,:,:,k));
end

C_tot_blk = zeros(N*L,N*L);

for l = 1:L
    C_tot_blk(1+(l-1)*N:l*N,1+(l-1)*N:l*N) = C_tot(:,:,l);
end

%Diagonal matrix with transmit powers and its square root
% Dp = diag(p);
% Dp12 = diag(sqrt(p));
% Dp112 = diag(sqrt(p1));

if ~strcmp(net_type, 'CENTRALIZED')
    V = nan(N,nbrOfRealizations,K,L);            
else
    V = nan(N*L,nbrOfRealizations,K);
end
W = V;
%% Go through all channel realizations
for n = 1:nbrOfRealizations
    if ~strcmp(net_type, 'CENTRALIZED') % CELLULAR or cell-free DISTRIBUTED
        %Go through all APs
        for l = 1:L
            
            %Extract channel estimate realizations from all UEs to AP l
            Hhatallj = reshape(Hhat(1+(l-1)*N:l*N,n,:),[N K]);
            
            V_MR = Hhatallj; %Compute MR combining
            V(:,n,:,l) = V_MR;
            if strcmp(combiner, 'M-MMSE') %Compute M-MMSE combining
                V(:,n,:,l) = p*((p*(Hhatallj*Hhatallj')+C_tot(:,:,l)+eyeN) \ V_MR); % A \ B = inv(A)*B,   A/B = A*inv(B)
            end
            W(:,n,:,l) = V(:,n,:,l) ./ vecnorm(V(:,n,:,l));
        end
        
    else %CENTRALIZED:  
        %Extract channel estimate realizations from all UEs to AP l
        Hhatallj = reshape(Hhat(:,n,:),[N*L K]);
        V_MR = Hhatallj; %Compute MR combining
        V(:,n,:) = V_MR;
        if strcmp(combiner, 'M-MMSE') %Compute M-MMSE combining
            V(:,n,:) = p*((p*(Hhatallj*Hhatallj')+C_tot_blk+eyeNL) \ V_MR); % A \ B = inv(A)*B,   A/B = A*inv(B)
        elseif strcmp(combiner, 'RZF') %Compute RZF combining
            V(:,n,:) = p*V_MR/(p*(V_MR'*V_MR)+eyeK);
        end
        W(:,n,:) = V(:,n,:) ./ vecnorm(V(:,n,:));
    end
    
end
