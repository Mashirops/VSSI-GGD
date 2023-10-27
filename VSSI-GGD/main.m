clc
clear
tic
%brainstorm %norgui
%%
algorithms = {'VSSI-GGD'};
WGN = 1; % Using white Gaussian Nosie (WGN = 1) or Human EEG noise (WGN  = 0);
LFNormlization = 0; % Whether normalize the LeadField Matrix
Uniform = 1; % Uniform/NonUniform Sources
VariousExtents = 0;
VariousSNRs = 1;
VariousSNIRs = 0;
VariousPatches = 0;
VariousCorrelation = 0;
VariousChannels = 0;
Test = 0;
if VariousExtents+VariousSNRs+VariousSNIRs+VariousPatches+VariousCorrelation+VariousChannels+Test ~= 1
    error('There will be one and only one scenario.');
end
% tic
%% Export the Channel,Cortex and HeadModel to workspace
if WGN
    channelselect=[1:32,34:42,44:64]; %Select EEG data
else
    channelselect=[1:32 34:42 44:59 61:63]; % for Real noise simulations
end

[sStudy, iStudy] = bst_get('StudyWithCondition','Subject01/Test');% import from brainstorm
index = 1;
bst_call(@export_matlab, {char(sStudy.Data(1,index).FileName)},'data');
%=======================Import the LFM and Cortex==================================%
[sSurface, iSurface] = bst_get('SurfaceFileByType',[],'Cortex');
% [sSurface, iSurface] = bst_get('SurfaceFileByType',2,'Cortex');
bst_call(@export_matlab, {char(sSurface.FileName)},'Cortex');
Atlas = Cortex.Atlas(2);

[sHeadModel] = bst_get('HeadModelForStudy', iStudy);
bst_call(@export_matlab, {char(sHeadModel.FileName)},'model');
Gain=model.Gain(channelselect,:);
GridLoc=model.GridLoc;
GridOrient=model.GridOrient;

Gain = bst_gain_orient(Gain,GridOrient);
clear GridOrient
[nSensor,nSource] = size(Gain);

%% Make output directory structure for imaging results (if doesn't already exist)
if VariousSNRs
   scenario = 'various SNRs';
   SNR1 = [0,3,5,10];
   SNIR1 = zeros(4,1);
   condition = SNR1';
   K = ones(4,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousSNIRs
   scenario = 'various SNIRs';
   SNR1 = zeros(4,1);
   SNIR1 = [0,3,5,10];
   condition = SNIR1';
   K = ones(4,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousExtents
   scenario = 'various extents';
   SNR1 = 5*ones(5,1);
   SNIR1 = 5*ones(5,1);
   condition = [1:5]';
   K = ones(5,1);
   DefinedArea = [2 5 10 18 32]'*1e-4*ones(1,max(K));%[0.5 4 8 14 22 32]'*1e-4*ones(1,2);% 38 48]'*1e-4;
elseif VariousChannels
   scenario = 'various channels';
   SNR1 = 5*ones(4,1);
   SNIR1 = 5*ones(4,1);
   condition = [62, 46, 32, 16]';
   K = ones(4,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousPatches
   scenario = 'various patches';
   SNR1 = 5*ones(4,1);
   SNIR1 = 5*ones(4,1);
   condition = [1:4]';
   K = [1,2,3,4];
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif VariousCorrelation
   scenario = 'various correlation coefficients';
   SNR1 = 5*ones(6,1);
   SNIR1 = zeros(6,1);
   condition = [0.1 0.3 0.5 0.7 0.9]';%(0:0.2:1)';
   K = 3*ones(6,1);
   DefinedArea = 8*1e-4*ones(size(condition,1),max(K));
elseif Test
    algorithms = {'VSSI-L2p'}
    ResultsLoading = 0;
    scenario = 'test';
    SNR1 = 10;
    SNIR1 = 10;
    condition = [1];
    K = 1;
    DefinedArea = 8*1e-4;
end

outpath = 'E:\result\VSSI-GGD\';
for i = 1 : size(condition,1)
    path{i} = fullfile(outpath,scenario,'\',num2str(condition(i)));
    if ~exist(path{i})
        mkdir(path{i});
    end
     if ~exist([path{i} '\' 'metrics.mat'], 'file')
         metrics = [];
         save([path{i} '\' 'metrics.mat'], 'metrics')
     end
end

%% Iteration
    dim = 0;
    Miter = 50 - dim;   
    Eccentricity = sqrt(sum(GridLoc.^2,2));
    Ec = find(Eccentricity > 70*1e-3);
    EstimatedArea_Mean = zeros(Miter,5);
for iter = 1:Miter  
    ind = randperm(numel(Ec));
    Thr = zeros(size(condition,1),numel(algorithms));
for iteration = 1:size(condition,1)
    fprintf('iter = %g, iteration = %g\n', iter,iteration)
    savepath = path{iteration};
    load ([path{iteration} '\' 'metrics'])
    SD = []; DLE = []; RMSE = []; AUC = []; PRE = []; REC = [];
    if any(ResultsLoading)
        SD = metrics.SD(iter,:); DLE = metrics.DLE(iter,:); RMSE = metrics.RMSE(iter,:); AUC = metrics.AUC(iter,:); PRE = metrics.PRE(iter,:); REC = metrics.REC(iter,:); nRMSE = metrics.nRMSE(iter,:); SE = metrics.SE(iter,:);
    end

%% Generate Simulated EEG Data
SNR = SNR1(iteration);
SNIR = SNIR1(iteration);
seedvox = Ec(ind(1:K(iteration)));
tau = [0.1 0.35 0.5 0.6];omega = [0.1 0.15 0.15 0.15];
f = [10 11 8 9];
Amp = 1e-8;
StimTime = find(abs(data.Time) == min(abs(data.Time)));
TimeLen = 300;
Time = data.Time(StimTime-0.5*TimeLen:StimTime+0.5*TimeLen-1);

OPTIONS.DefinedArea    = DefinedArea(iteration,:);
OPTIONS.seedvox        = seedvox;
OPTIONS.frequency      = f;
OPTIONS.tau            = tau;
OPTIONS.omega          = omega;
OPTIONS.Amp            = Amp;
OPTIONS.GridLoc        = GridLoc;
if VariousPatches
OPTIONS.MixMatrix = [1    0    0     0;
                     0    1    0     0;
                     0    0    1     0;
                     0    0    0     1;
                     0    0    .5     .5];
elseif VariousCorrelation
    xi = condition(iteration);
OPTIONS.MixMatrix = [1            0            0          0;
                     xi     sqrt(1-xi^2)       0          0;
                     xi           0        sqrt(1-xi^2)   0;
                     0            0            1          0];
end
OPTIONS.uniform       = Uniform;
OPTIONS.WGN           = WGN;
OPTIONS.SNR           = SNR;
OPTIONS.SNIR          = SNIR;
OPTIONS.ar            = 0;
OPTIONS.params(:,:,1) = [ 0.8    0    0 ;
                            0  0.9  0.5 ;
                          0.4    0  0.5];

OPTIONS.params(:,:,2) = [-0.5    0    0 ;
                            0 -0.8    0 ;
                            0    0 -0.2];

OPTIONS.noisecov      = [ 0.3    0    0 ;
                            0    1    0 ;
                            0    0  0.2];

if ~any(ResultsLoading)
    [Data,s_real,Result] = Simulation_Data_Generate (Gain,Cortex,Time,OPTIONS);
    ActiveVoxSeed = Result.ActiveVoxSeed;
else
    % %================= Recording from previous Results ====================%
    load ([savepath '\' 'result' num2str(iter+dim)]);
    Data = Result.B;
    s_real = Result.real;
    seedvox = Result.seedvox;
    Result.B = Data;
    [~, VertArea] = tess_area(Cortex.Vertices, Cortex.Faces);
    AreaDef = DefinedArea(iteration,:);
    ActiveVox = [];
    for k = 1:numel(seedvox)
        ActiveVoxSeed{k} = PatchGenerate(seedvox(k),Cortex.VertConn,VertArea,AreaDef(k));
        ActiveVox = union(ActiveVoxSeed{k},ActiveVox);
    end
    StimTime = size(Data,2)/2 + 1;
end
% %=======================================================================%
fprintf('Actual SNR is %g\n',20*log10(norm(Gain*s_real,'fro')/norm(Data-Gain*s_real,'fro')));
corr(s_real(seedvox,StimTime+1:end)','type','Pearson')
 
%======================================%
%         Leadfield Matrix normalization
%=====================================%
if LFNormlization
    LfvW = sqrt(sum(Gain_scale.^2,1).^0.3);
    Gain_scale = Gain.*kron(ones(nSensor,1),1./LfvW);
end
%% Whiten measurements and lead field matrix
Bpre = B(:,1:StimTime);
Cov_n = NoiseEstimate(B,StimTime);
NoiseMethod = 'median';%'reg';%'none';%{'shrink', 'reg', 'diag', 'none', 'median'};
FourthMoment = Bpre*Bpre'./StimTime;nSamples = StimTime;
NoiseReg = 0.1;
[Cov,W] = truncate_and_regularize_covariance(Cov_n,NoiseMethod,'EEG',NoiseReg,FourthMoment,nSamples);
L = W*Gain_scale;
B = W*B;
if ~any(ResultsLoading)
    Result.Whiter = W;
end
clear W;
%%  SVD for TBFs
[Dic] = TBFSelection(B,0);%,'threshold','Permutation');'Kaiser'
%% Source Estimation
MetricsInterval = [];

%% Channel select
if VariousChannels
    if ~any(ResultsLoading)
        cnumber = condition(iteration);
        load(['C:\Users\guest0\Desktop\GGD-fig\channel\',num2str(cnumber),'c.mat']);
        channels(channels>33 & channels<43) = channels(channels>33 & channels<43) - 1;
        channels(channels>43) = channels(channels>43) - 2;
        B = B(channels,:);
        L = L(channels,:);
        Result.channels = channels;
    else
        channels = Result.channels;
        B = B(channels,:);
        L = L(channels,:);

    end
end

%% VSSI-GGD-ADMM
Weight = logspace(-4,0,20);
variation = 'Sparse+Variation';
opts.sparse = 0.67;%Weight(14);
if any(strcmpi('VSSI-GGD', algorithms))
    MethodIndex = find(strcmpi('VSSI-GGD', algorithms)~=0);
    Edge = VariationEdge(Cortex.VertConn);
    if ~ResultsLoading(MethodIndex)
        [S_GGD] = VSSI_GGD_ADMM(B*Dic',L,Cortex.VertConn,opts.sparse,'transform',variation,'p',0.8,'tol',1e-6,'roupar',1e17);
        S_vssiggd = S_GGD{end}*Dic;
        if LFNormlization
            S_vssiggd = repmat(1./LfvW',1,size(S_vssiggd,2)).*S_vssiggd;
        end
        S_vssiggd = S_vssiggd*ratio;
    else
        S_vssiggd = Result.VSSIGGD;
    end
    
    [SD(1,MethodIndex),DLE(1,MethodIndex),RMSE(1,MethodIndex),nRMSE(1,MethodIndex),SE(1,MethodIndex),PRE(1,MethodIndex),REC(1,MethodIndex)]...
        = PerformanceMetric(GridLoc,S_vssiggd(:,StimTime:end),s_real(:,StimTime:end),ActiveVoxSeed);%,'interval',MetricsInterval);
    Roc = ROCextent(s_real(:,StimTime:end),S_vssiggd(:,StimTime:end),Cortex,seedvox);
    AUC(1,MethodIndex) = median(Roc.mean);
    EV(1,MethodIndex) =  1 - norm(B*ratio - L*S_vssiggd,'fro')^2/norm(B*ratio,'fro')^2;
    Result.VSSIGGD = S_vssiggd;
end 

%% Save Results on Disk
method = 1:numel(algorithms);
metrics.AUC(dim+iter,method) = AUC;%(method);
metrics.SD(dim+iter,method) = SD;%(method); 
metrics.DLE(dim+iter,method) = DLE;%(method);
metrics.RMSE(dim+iter,method) = RMSE;%(method);
metrics.nRMSE(dim+iter,method) = nRMSE;%(method);
metrics.SE(dim+iter,method) = SE;%(method);
metrics.PRE(dim+iter,method) = PRE;
metrics.REC(dim+iter,method) = REC;


save_file_name=[savepath '\' 'result' num2str(dim+iter)];
save (save_file_name,'Result');
save([savepath '\' 'metrics.mat'], 'metrics')


% toc
end
end
save([savepath '\' 'GridLoc.mat'], 'GridLoc')
save([savepath '\' 'Cortex.mat'], 'Cortex')
save([savepath '\' 'methods.mat'], 'algorithms')
runningtime = toc;