function Thr = ThresholdSelect(S,varargin)
% Discription: Automated choose the threshold for E/MEG source imaging
% results.
% Reference: Otsu N. A threshold selection method from gray-level histograms[J]. Automatica, 1975, 11(285-296): 23-27.
% By Liu Ke on 2015/10/13
ifplot = 0;
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'ifplot'   
                ifplot = varargin{i+1}; 
        end
    end
end
S = sqrt(sum(S.^2,2));
S = abs(S)./max(abs(S));
[a,b] = hist(S,100);
Cost = zeros(numel(a),1);
P = a./numel(S);
mu_T = sum(P.*b);
for k = 1: numel(a)
    w_k = sum(P(1:k));
    mu_k = sum(P(1:k).*b(1:k));
    Cost(k) = (mu_T*w_k - mu_k)^2/w_k/(1-w_k);
end
sigma_T = sum((b -mu_T).^2.*P);
yita = Cost./sigma_T;
[w,idx] = max(yita);
if ifplot
%     figure
%     plot(b,Cost);
    figure
    plot(b,yita);
end
Thr = b(idx);