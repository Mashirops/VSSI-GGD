function Sigma = NoiseEstimate(B,StimTime)
% Discription: Automatedly estimate the covariance of E/MEG data using
% cross validation
% Reference: Engemann D A, Gramfort A. Automated model selection in 
%            covariance estimation and spatial whitening of MEG and EEG signals[J]. NeuroImage, 2015, 108: 328-342.
 K = 3; % Number of folds in cross..validation
 if isempty(StimTime)
     Bpre = B;
     Bpost = [];
 else
     Bpre = B(:,1:StimTime);
     Bpost = B(:,StimTime+1:size(B,2));
 end
 
 [indices] = crossvalind('Kfold',size(Bpre,2),K); % Find split

Likelihood = @(B,Sigma) -trace(B*B'/Sigma)/(2*size(B,2)) - (log(det(10^(-floor(log10(mean(diag(Sigma)))))*Sigma)) + numel(Sigma(:))*floor(log10(min(diag(Sigma))))*log(10))/2;
alpha_all = linspace(0.01,1,30);
L = zeros(1,numel(alpha_all));
Alpha = zeros(K,1);
for iter = 1:K
    B_train = Bpre(:,indices==iter);
    B_test = [Bpre(:,indices~=iter),Bpost];
    
   C = B_train*B_train'./size(B_train,2);
   mu = mean(diag(C));
    
    for i = 1:numel(alpha_all)
        Sigma = (1 - alpha_all(i))*C + alpha_all(i)*mu*eye(size(B,1));
        L(i) = Likelihood(B_test,Sigma);
    end
    [~, imax] = max(L);
    Alpha(iter) = alpha_all(imax);
end
Alpha = mean(Alpha);
C = Bpre*Bpre'./size(Bpre,2);
mu = mean(diag(C));
Sigma = (1 - Alpha)*C + Alpha*mu*eye(size(B,1));