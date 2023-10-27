%% =========== Covariance Truncation and Regularization
function [Cov,iW] = truncate_and_regularize_covariance(Cov,Method,Type,NoiseReg,FourthMoment,nSamples)
%         |- NoiseCovMat        : Noise covariance structure
%         |   |- NoiseCov       : Noise covariance matrix
%         |   |- FourthMoment   : Fourth moment (F^2 * F^2'/n)
%         |   |- nSamples       : Number of times samples used to compute those measures

% Cov is the covariance matrix, to be regularized using Method
% Type is the sensor type for display purposes
% NoiseReg is the regularization fraction, if Method "reg" selected
% FourthMoment and nSamples are used if Method "shrinkage" selected

VERBOSE = true; % be talkative about what's happening

% Ensure symmetry
Cov = (Cov + Cov')/2;

% Note,impossible to be complex by above symmetry check
% Decompose just this covariance.
[Un,Sn2] = svd(Cov,'econ');
Sn = sqrt(diag(Sn2)); % singular values
tol = length(Sn) * eps(single(Sn(1))); % single precision tolerance
Rank_Noise = sum(Sn > tol);

if VERBOSE,
    fprintf('BST_INVERSE > Rank of the %s channels, keeping %.0f noise eigenvalues out of %.0f original set\n',...
        Type,Rank_Noise,length(Sn));
end

Un = Un(:,1:Rank_Noise);
Sn = Sn(1:Rank_Noise);

% now rebuild the noise covariance matrix with just the non-zero components
Cov = Un*diag(Sn.^2)*Un'; % possibly deficient matrix now

% With this modality truncated, see if we need any additional
% regularizations, and build the inverse whitener

if VERBOSE
    fprintf('BST_INVERSE > Using the ''%s'' method of covariance regularization.\n',Method);
end

switch(Method) % {'shrink', 'reg', 'diag', 'none', 'median'}
    
    case 'none'
        %  "none" in Regularization means no  regularization was applied to the computed Noise Covariance
        % Matrix. Do Nothing to Cw_noise
        iW = Un*diag(1./Sn)*Un'; % inverse whitener
        if VERBOSE,
            fprintf('BST_INVERSE > No regularization applied to covariance matrix.\n');
        end
        
        
    case 'median'
        if VERBOSE,
            fprintf('BST_INVERSE > Covariance regularized by flattening tail of eigenvalues spectrum to the median value of %.1e\n',median(Sn));
        end
        Sn = max(Sn,median(Sn)); % removes deficient small values
        Cov = Un*diag(Sn.^2)*Un'; % rebuild again.
        iW = Un*diag(1./Sn)*Un'; % inverse whitener
        
    case 'diag'
        Cov = diag(diag(Cov)); % strip to diagonal
        iW = diag(1./sqrt(diag(Cov))); % inverse of diagonal
        if VERBOSE,
            fprintf('BST_INVERSE > Covariance matrix reduced to diagonal.\n');
        end
        
    case 'reg'
        % The unit of "Regularize Noise Covariance" is as a percentage of
        % the mean variance of the modality.
        
        % Ridge Regression:
        RidgeFactor = Sn2(1) * NoiseReg ; % percentage of max
%         RidgeFactor = mean(diag(Sn2)) * NoiseReg ; % percentage of mean
        Cov = Cov + RidgeFactor * eye(size(Cov,1));
        iW = Un*diag(1./(Sn + sqrt(RidgeFactor)))*Un'; % inverse whitener
        
        if VERBOSE,
            fprintf('BST_INVERSE > Diagonal of %.1f%% of largest eigenvalue added to covariance matrix.\n',NoiseReg * 100);
        end
        
        
    case 'shrink'
        % Method of Ledoit, recommended by Alexandre Gramfort
        
        % Need to scale the Fourth Moment for the modalities
       
        % use modified version of cov1para attached to this function
        % TODO, can we adjust this routine to handle different numbers of
        % samples in the generation of the fourth order moments
        % calculation? As of August 2016, still relying on a single scalar
        % number.
        [Cov,shrinkage]=cov1para_local(Cov,FourthMoment,nSamples);
        if VERBOSE,
            fprintf('\nShrinkage factor is %f\n\n',shrinkage)
        end
        % we now have the "shrunk" whitened noise covariance
        % Recalculate
        [Un,Sn2] = svd(Cov,'econ');
        Sn = sqrt(diag(Sn2)); % singular values
        tol = length(Sn) * eps(single(Sn(1))); % single precision tolerance
        Rank_Noise = sum(Sn > tol);
        
        if VERBOSE,
            fprintf('BST_INVERSE > Ledoit covariance regularization, after shrinkage, rank of the %s channels, keeping %.0f noise eigenvalues out of %.0f original set\n',...
                Type,Rank_Noise,length(Sn));
        end
        
        Un = Un(:,1:Rank_Noise);
        Sn = Sn(1:Rank_Noise);
        
        % now rebuild the noise covariance matrix with just the non-zero
        % components
        Cov = Un*diag(Sn.^2)*Un'; % possibly deficient matrix now
        
        iW = Un*diag(1./Sn)*Un'; % inverse whitener
        
        
    otherwise
        error(['BST_INVERSE > Unknown covariance regularization method: NoiseMethod="' Method '"']);
        
end % method of regularization


% Note the design of full rotating whiteners. We don't expect dramatic reductions in
% rank here, and it's convenient to rotate back to the original space.
% Note that these whitener matrices may not be of full rank.

end
%% ======== Ledoit Shrinkage
% Modified to use the precalculated stats

function [sNoiseCov,shrinkage]=cov1para_local(NoiseCov,FourthMoment,nSamples)
% Based on Ledoit's "cov1para" with some modifications
%   x is t x n, returns
%   sigma n x n
%
% shrinkage is the final computed shrinkage factor, used to weight the
%  i.i.d. prior vs the sample estimate. If shrinkage is specified, it is
%  used on input; else, it's computed.

% Original code from
% http://www.ledoit.net/cov1para.m
% Original Ledoit comments:
% function sigma=cov1para(x)
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards one-parameter matrix:
%    all variances are the same
%    all covariances are zero
% if shrink is specified, then this const. is used for shrinkage

% Based on
% http://www.ledoit.net/ole1_abstract.htm
% http://www.ledoit.net/ole1a.pdf (PDF of paper)
%
% A Well-Conditioned Estimator for Large-Dimensional Covariance Matrices
% Olivier Ledoit and Michael Wolf
% Journal of Multivariate Analysis, Volume 88, Issue 2, February 2004, pages 365-411
%
% Abstract
% Many economic problems require a covariance matrix estimator that is not
% only invertible, but also well-conditioned (that is, inverting it does
% not amplify estimation error). For large-dimensional covariance matrices,
% the usual estimator - the sample covariance matrix - is typically not
% well-conditioned and may not even be invertible. This paper introduces an
% estimator that is both well-conditioned and more accurate than the sample
% covariance matrix asymptotically. This estimator is distribution-free and
% has a simple explicit formula that is easy to compute and interpret. It
% is the asymptotically optimal convex combination of the sample covariance
% matrix with the identity matrix. Optimality is meant with respect to a
% quadratic loss function, asymptotically as the number of observations and
% the number of variables go to infinity together. Extensive Monte-Carlo
% confirm that the asymptotic results tend to hold well in finite sample.


% Original Code Header, updated to be now from (2014)
% http://www.econ.uzh.ch/faculty/wolf/publications/cov1para.m.zip
%
% x (t*n): t iid observations on n random variables
% sigma (n*n): invertible covariance matrix estimator
%
% Shrinks towards constant correlation matrix
% if shrink is specified, then this constant is used for shrinkage
%
% The notation follows Ledoit and Wolf (2003, 2004)
% This version 04/2014
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file is released under the BSD 2-clause license.
%
% Copyright (c) 2014, Olivier Ledoit and Michael Wolf
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
% 1. Redistributions of source code must retain the above copyright notice,
% this list of conditions and the following disclaimer.
%
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
% IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
% THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
% PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
% CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
% EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
% PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
% PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
% LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
% NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
% SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Wolf's site now,
% http://www.econ.uzh.ch/faculty/wolf/publications/cov1para.m.zip
% some differences from original Ledoit that confirm Mosher's
% original re-coding.

% % de-mean returns
% [t,n]=size(x);
% meanx=mean(x);
% x=x-meanx(ones(t,1),:);

% compute sample covariance matrix
% Provided
% NoiseCov=(1/t).*(x'*x);

% compute prior
n=size(NoiseCov,1); % number of channels
meanvar=mean(diag(NoiseCov));
prior=meanvar*eye(n); % Note, should be near identity by our pre-whitening

% what we call p
%y=x.^2;
%phiMat=y'*y/t - NoiseCov.^2;
phiMat = FourthMoment - NoiseCov.^2;
phi=sum(sum(phiMat));

% what we call r is not needed for this shrinkage target

% what we call c
gamma=norm(NoiseCov-prior,'fro')^2;

% compute shrinkage constant
kappa=phi/gamma;
% ensure bounded between zero and one
shrinkage=max(0,min(1,kappa/nSamples));

% compute shrinkage estimator
sNoiseCov=shrinkage*prior+(1-shrinkage)*NoiseCov;

end