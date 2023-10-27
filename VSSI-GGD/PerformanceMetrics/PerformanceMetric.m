% Discription: This script is to calculate the perfomance metric for a
% specified localization method
% Two metrics:
% (1)spatial dispersion (SD), which measures how source locations are spatially
%   blurred;
% (2) distance of localization error (DLE), which measures the 
%   averaged error distance between the true source
%   locations and the maximum of source estimates.
% (3) RMSE: To assess the ability to recover the theoretical distribution of current s_real with an accurate amplitude
% (4) SE: Shape error, the error between the normalized simulated and estimated sources
% Input : GridLoc : 3D locacation of each dipole within the brain
%         s: reconstructed currents
%         s_real: simulated currents(Ground Truth)
% Output: SD and DLE

% Reference : (1)Chang W, Nummenmaa A, Hsieh J, Lin F. 
% Spatially sparse source cluster modeling by compressive neuromagnetic tomography. 
% Neuroimage 53(1): 146--160, 2010.
% (2) Grova C, Daunizeau J, Lina JM, Benar CG, Benali H, Gotman J. 
% Evaluation of EEG localization methods using realistic simulations of interictal  spikes. 
% Neuroimage 29(3): 734-53, 2006.

% Author : Ke Liu
% Date : 2013/11/15


function [SD,DLE,RMSE,SE,SC,PRE,REC] = PerformanceMetric(GridLoc,s,s_real,ActiveVoxSeed,varargin)

% SD = 0; DLE = 0;
TimeInterval = 1:size(s,2);
for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'interval'   
                TimeInterval = varargin{i+1}; 
            otherwise
                error(['Unrecognized parameter: ''' varargin{i} '''']);
        end
end
if isempty(TimeInterval)
    TimeInterval = 1:size(s,2);
end
          keeplist = (1:size(s,1))';          
          P = sqrt(sum(s(:,TimeInterval).^2,2));
%           index = (P~=0);
          index = (P>=0.1*max(P));
          keeplist = keeplist(index);
          ActiveVox = find(mean(s_real(:,TimeInterval).^2,2)~=0);
          I = SpatialNeigh(ActiveVox,GridLoc,keeplist);
          SD = SpatialDispersion(s(:,TimeInterval),I,ActiveVox,GridLoc);
          DLE = DisLE(s(:,TimeInterval),I,ActiveVox,GridLoc);

          
%           s = s.*norm(s_real,'fro')/norm(s,'fro'); % Normalize

%           s = s./norm(s,'fro');
%           s_real = s_real./norm(s_real,'fro');
% s = s./max(abs(s(:)));
% s_real = s_real./max(abs(s_real(:)));
          %% MSE within the whole cortex
          RMSE = norm(s-s_real,'fro')^2/norm(s_real,'fro')^2;
          %% MSE within the simulated Patch 
          %Shape error of the sum of all simulated patches
%           sn = s./repmat( ( max( abs(s') ) )' ,1,size(s,2));
%           s_realn = s_real./repmat(( max( abs(s_real') ) )',1,size(s_real,2));
%           SE = 0;
%           for i = 1:numel(ActiveVoxSeed)
%                temp = mean(sn(ActiveVoxSeed{i},:)); temp = temp./max(abs(temp));
%                temp1 = mean(s_realn(ActiveVoxSeed{i},:)); temp1 = temp1./max(abs(temp1));
%                SE = sqrt(sum(temp - temp1).^2/size(s,2)) + SE;
%           end
          
          
          s = s./norm(s,'fro');
          s_real = s_real./norm(s_real,'fro');
          SE = norm(s-s_real,'fro')^2/norm(s_real,'fro')^2;
          %SE = norm(s(ActiveVox,:) - s_real(ActiveVox,:),'fro')^2/norm(s_real(ActiveVox,:),'fro')^2;
          %% Correlation coefficient between mean simulated and estimated time courses within the active areas
          c = zeros(numel(ActiveVoxSeed),1);
          for i = 1:numel(ActiveVoxSeed)
              c(i) = corr2(mean(s(ActiveVoxSeed{i},:),1),s_real(ActiveVoxSeed{i}(1),:));
          end
          SC = mean(c);
%           smean  = mean(s(ActiveVox,:))./max(mean(s(ActiveVox,:)));
%           srealmean = mean(s_real(ActiveVox,:))./max(mean(s_real(ActiveVox,:)));
%           SE = corr2(smean,srealmean);%norm(smean - srealmean,'fro')^2;
 %% Precision and Recall
          T = ThresholdSelect(s);
          GT = sqrt(sum(s_real.^2,2)); GT(GT~=0) = 1;
          RC = sqrt(sum(s.^2,2));
          RC = abs(RC)./max(abs(RC)); 
          RC(RC<T) = 0; RC(RC>=T) = 1;
          OL = GT.*RC;
          PRE = sum(OL)/sum(RC);
          REC = sum(OL)/sum(GT);
end


        
function I = SpatialNeigh(ActiveVox,GridLoc,keeplist)  
% only the dipoles whose power are larger than 1% of the maximum power is
% taken into consideration;
% keeplist: the index number of dipoles whose power is larger than 10% of
% the maximum power
    I = cell(numel(ActiveVox),1);      
    for j =1:numel(keeplist)
%          distance = zeros(numel(ActiveVox),1);
%        for k =1 :numel(ActiveVox)  
%            distance(k) = norm(GridLoc(keeplist(j),:)-GridLoc(ActiveVox(k),:));
%        end
       distance = sqrt(sum((GridLoc(ActiveVox,:)-repmat(GridLoc(keeplist(j),:),numel(ActiveVox),1)).^2,2));
       
       [~,IX] = sort(distance,'ascend');
       I{IX(1)} = union(I{IX(1)},keeplist(j));
%           end
    end
end   
    
    
    
function SD = SpatialDispersion(s,I,ActiveVox,GridLoc) 
     SD = 0;
    for i = 1:numel(ActiveVox)
        if numel(I{i})~=0
            for j = 1:numel(I{i})
                SD = SD + norm(GridLoc(I{i}(j),:)-GridLoc(ActiveVox(i),:))^2*sum(s(I{i}(j),:).^2);
            end
        end
    end
    SD = SD/sum(sum(s.^2));
    SD = sqrt(SD);
end
  
    
    
 function DLE = DisLE(s,I,ActiveVox,GridLoc)    
    J = [];
    for i =1 :numel(ActiveVox)
        if numel(I{i})~=0
            J = union(J,i);
        end
    end
    distance = zeros(numel(J),1);
    for i = 1:numel(J)
            power = sum(s(I{J(i)},:).^2,2);
            [~,IX] = sort(power,'descend');
            distance(i) = norm(GridLoc(ActiveVox(J(i)),:)-GridLoc(I{J(i)}(IX(1)),:));
    end
    DLE = mean(distance);
 end
