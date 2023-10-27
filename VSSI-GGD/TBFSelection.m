function [TBFs] = TBFSelection(B,ifplot,varargin)

ThresholdMethod = 'Kaiser';
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'threshold'
                ThresholdMethod = varargin{i+1}; 
        end
    end
end
%% for permutation
if strcmpi(ThresholdMethod, 'Permutation')
    nperms = 1000;
    permeigvals = zeros(nperms,size(B,1));
    permsigvals = zeros(nperms,size(B,1));
    for permi=1:nperms
        % random data by randomizing ERP time points/channels
        randdat = reshape(B(randperm(numel(B ))),size(B));
        [junk,v,junk] = svd(randdat,'econ');
        
        eigvals = diag(v).^2;
        permeigvals(permi,:) = eigvals./sum(eigvals);
        sigvals = diag(v);
        permsigvals(permi,:) = sigvals./sum(sigvals);
    end
    
    for i = 1:size(B,1)
        [a,junk] = sort(permeigvals(:,i),'descend');
        threshold(i) = a(nperms*0.01);
    end
end
%% Kaiser criterion
if strcmpi(ThresholdMethod, 'Kaiser')
    threshold = 1/size(B,1);
end
%% for original data
[junk,v,d] = svd(B,'econ');
Eigvals = diag(v).^2;
NEigvals = Eigvals./sum(Eigvals);
%% TBFselect
KSVD = 6;%sum(NEigvals' >= threshold);
TBFs = d(:,1:KSVD)';
%% Figure Plot
if ifplot == 1
    if strcmpi(ThresholdMethod, 'Permutation')
        figure
        plot(NEigvals(1:min(30,size(B,1)) ),'-o','linewidth',1);
        hold on
        plot(threshold(1:min(30,size(B,1)) ),'r-p','markerface','w','linewidth',1)
        legend('Normalized Eignvalues','Threshold(permutation)')
        hold on
        plot((KSVD + .5)*ones(1,2),[0,max(NEigvals)],'k-.','linewidth',1)
        text(KSVD + .5,0.7*max(NEigvals),['\leftarrow K= ',num2str(KSVD)])
        xlabel('Componenet number');ylabel('Normalized Eignvalues')
    end
    
    
    if strcmpi(ThresholdMethod, 'Kaiser')
        figure
        plot(NEigvals(1:min(30,size(B,1)) ),'-o','linewidth',1);
        hold on
        plot([0,min(30,size(B,1))],1/size(B,1)*ones(1,2),'g-.','linewidth',2)
        legend('Normalized Eignvalues','Threshold (analytic)')
        hold on
        plot((KSVD + 0.5)*ones(1,2),[0,max(NEigvals)],'k-.','linewidth',1)
        text(KSVD + 0.5,0.7*max(NEigvals),['\leftarrow K= ',num2str(KSVD)])
        xlabel('Componenet number');ylabel('Normalized Eignvalues')
    end
    
end