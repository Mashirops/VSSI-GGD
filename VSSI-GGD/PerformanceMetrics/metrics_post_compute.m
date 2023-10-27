 clc
 clear
 
 path{1} = 'E:\result\STBF\Simulations1';
 
% path{1} = 'E:\result\GLM\various extents\1';
% path{2} = 'E:\result\GLM\various extents\2';
% path{3} = 'E:\result\GLM\various extents\3';
% path{4} = 'E:\result\GLM\various extents\4';
% path{5} = 'E:\result\GLM\various extents\5';

% path{1} = 'E:\result\GLM\various patches\1';
% path{2} = 'E:\result\GLM\various patches\2';
% path{3} = 'E:\result\GLM\various patches\3';
% path{4} = 'E:\result\GLM\various patches\4';

% path{1} = 'E:\result\GLM\various SNRs\1';
% path{2} = 'E:\result\GLM\various SNRs\2';
% path{3} = 'E:\result\GLM\various SNRs\3';
% path{4} = 'E:\result\GLM\various SNRs\4';

% path{1} = 'E:\result\TBMRF\various patches\1';
% path{2} = 'E:\result\TBMRF\various patches\2';
% path{3} = 'E:\result\TBMRF\various patches\3';
% path{4} = 'E:\result\TBMRF\various patches\4';
% load (['E:\result\TBMRFInf\various extents' '\' 'Gain' ])
for Iteration = 1:numel(path)
    savepath = path{Iteration};
load ([savepath '\' 'GridLoc']);
load ([savepath '\' 'Cortex']);
% load ([savepath '\' 'real']);

dim = 0;
Miter = 50 - dim;
SD = zeros(Miter,5);DLE = SD; RMSE = SD; %AUC = SD; OSNIR = SD;
for iter = 1:Miter
 disp(['Iteration:   ',num2str(iter)])  
load ([savepath '\' 'result' num2str(iter+dim)]);
s_real = Result.real;
seedvox = Result.seed;%vox;
% seedvox = 10833;



    B = Result.B;
    W =  Result.Whiter;
%     StimTime = 251;
%     Bpre = B(:,1:StimTime);
%     Cov_n = (Bpre - repmat(mean(Bpre,2),1,StimTime)) * (Bpre - repmat(mean(Bpre,2),1,StimTime))'/(StimTime - 1);
%     [V,D] = eig(Cov_n); 
%     D = diag(D); 
%     [D,I] = sort(D,'descend'); 
%     V = V(:,I);
%     D = 1./D;
%     W = diag(sqrt(D)) * V';
%     clear V D I


for i = 1:numel(seedvox)
    ActiveVox{i} = seedvox;
end
[SD(iter,1),DLE(iter,1),RMSE(iter,1),nRMSE(iter,1),SE(iter,1)] ...
    = PerformanceMetric(GridLoc,Result.STBF,s_real,ActiveVox);

[SD(iter,2),DLE(iter,2),RMSE(iter,2),nRMSE(iter,2),SE(iter,2)] ...
    = PerformanceMetric(GridLoc,Result.SBFKernel*W*B,s_real,ActiveVox);

[SD(iter,3),DLE(iter,3),RMSE(iter,3),nRMSE(iter,3),SE(iter,3)] ...
    = PerformanceMetric(GridLoc,Result.sLORETA*W*B,s_real,ActiveVox);

[SD(iter,4),DLE(iter,4),RMSE(iter,4),nRMSE(iter,4),SE(iter,4)] ...
    = PerformanceMetric(GridLoc,Result.LORETA*W*B,s_real,ActiveVox);

% [SD(iter,1),DLE(iter,1),RMSE(iter,1)] = PerformanceMetric(GridLoc,Result.Mean,s_real);
%  A = ROCextent(s_real,Result.TBMRF,Cortex,seedvox);
%  AUC(iter,1) = median(A.mean);


% [SD(iter,2),DLE(iter,2),RMSE(iter,2)] = PerformanceMetric(GridLoc,Result.Map,s_real);
% A = ROCextent(s_real,Result.TNAdaptive,Cortex,seedvox);
% AUC(iter,1) = median(A.mean);
 
% [SD(iter,3),DLE(iter,3),RMSE(iter,3)] = PerformanceMetric(GridLoc,Result.wMNE*W*B,s_real); 
% A = ROCextent(s_real,Result.wMNE*W*B,Cortex,seedvox);
% AUC(iter,2) = median(A.mean);

% [SD(iter,4),DLE(iter,4),RMSE(iter,4)] = PerformanceMetric(GridLoc,Result.LORETA*W*B,s_real);
% A = ROCextent(s_real,Result.LORETA*W*B,Cortex,seedvox);
% AUC(iter,3) = median(A.mean);

% [SD(iter,4),DLE(iter,4),RMSE(iter,4)] = PerformanceMetric(GridLoc,Result.SBLkernel*W*Result.B,s_real);
% A = ROCextent(s_real,Result.SBLkernel*W*B,Cortex,seedvox,'flag',2);
% AUC(iter,4) = median(A.mean);

% [SD(iter,5),DLE(iter,5),RMSE(iter,5)] = PerformanceMetric(GridLoc,Result.vbglm,s_real);
%  A = ROCextent(s_real,Result.vbglm,Cortex,seedvox);
%  AUC(iter,5) = median(A.mean);

% [SD(iter,5),DLE(iter,5),RMSE(iter,5)] = PerformanceMetric(GridLoc,Result.Sparse,s_real);
% Roc = ROCextent(s_real,Result.Sparse,Cortex,seedvox,'flag',2);
% AUC(iter,5) = median(Roc.mean);

% OSNIR(iter,1) = 10*log10(norm(Gain*s_real,'fro')^2/norm(Gain*Result.TBMRFInf-Gain*s_real,'fro')^2);
% OSNIR(iter,2) = 10*log10(norm(Gain*s_real,'fro')^2/norm(Gain*Result.TBMRF-Gain*s_real,'fro')^2);
% OSNIR(iter,3) = 10*log10(norm(Gain*s_real,'fro')^2/norm(Gain*Result.wMNE*W*B-Gain*s_real,'fro')^2);
% OSNIR(iter,4) = 10*log10(norm(Gain*s_real,'fro')^2/norm(Gain*Result.LORETA*W*B-Gain*s_real,'fro')^2);
% OSNIR(iter,5) = 10*log10(norm(Gain*s_real,'fro')^2/norm(Gain*Result.SBLkernel*W*Result.B-Gain*s_real,'fro')^2);
end

   load ([savepath '\' 'metrics'])

  method = [1:4];%[3:5];
  metrics.SD(1+dim:Miter+dim,method) = SD(:,method);
  metrics.DLE(1+dim:Miter+dim,method) = DLE(:,method);
  metrics.nRMSE(1+dim:Miter+dim,method) = RMSE(:,method); 
%   metrics.AUC(1+dim:Miter+dim,method) = AUC(:,method);
% metrics.OSNIR(1+dim:Miter+dim,method) = OSNIR(:,method);
  save([savepath '\' 'metrics.mat'], 'metrics')
end