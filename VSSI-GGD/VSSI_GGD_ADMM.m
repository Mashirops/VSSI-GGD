function [S,cost] = VSSI_GGD_ADMM(B,L,VertConn,sparseweight,varargin)
%% Description: Reconstructed extended sources based on Generalized Gaussian prior
% Data Model:
% B = LS + epsilon;
% p(epsilon(:,t)) = N(0,I);
% U = V*S; variation sources

% Input:
%         B(d_b x T):               M/EEG Measurement
%         L(d_b x d_s):             Leadfiled Matrix
%         VertConn:                 Cortex Connectivity Condition
%         sparseweight:             sparseweight for the sparse variation 稀疏变化的稀疏权值
%                                   (typically = 0.01)

% Output:
%         S:                        Estimated Sources

%% 
[Nsensor,Nsource] = size(L);
T = size(B,2);                  %T sample points
V = VariationEdge(VertConn);
tol = 1e-3;                         %tolerance
QUIET = 1;
rou_update = 1;
cost = 0;
rou = 1e15;
p = 0.7;

% get input argument values
if(mod(length(varargin),2)==1)
    error('Optional parameters should always go by pairs\n');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'transform'
                transform = varargin{i+1};
            case 'tol'
                tol = varargin{i+1};
            case 'roupar'
                rou = varargin{i+1};
            case 'p'
                p = varargin{i+1};
        end
    end
end       
  Edge = VariationEdge(VertConn);      
if strcmp(transform, 'Variation')
    V = VariationEdge(VertConn);
elseif strcmp(transform, 'Laplacian')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
elseif strcmp(transform, 'Laplacian+Variation')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
    V = [opts.laplacian*V;VariationEdge(VertConn)];
elseif strcmp(transform,'Sparse+Laplacian')
    NVertConn = sum(VertConn,2);
    V = bsxfun(@minus,spdiags(ones(Nsource,1),0,Nsource,Nsource),bsxfun(@times,bsxfun(@rdivide,VertConn,NVertConn),0.95*ones(Nsource,1)));
    V = [sparseweight*sparse(1:Nsource,1:Nsource,1);V];
elseif strcmp(transform,'Sparse+Variation')
    V = [sparseweight*sparse(1:Nsource,1:Nsource,1);VariationEdge(VertConn)];
elseif strcmp(transform, 'Sparse')
    V = sparse(1:Nsource,1:Nsource,1);
end




Max_iter = 10;
ADMM_iter = 400;
U_iter = 10;

% Initial values
TMNE = MNE(B,[],L,[],'MNE');
S_MNE = TMNE*B;

S = S_MNE;%zeros(Nsource,T);%
Z = zeros(size(V,1),T);
U = V*S;%zeros(size(V,1),T);    %V*S;% 
U_old = U;

tmp = sqrt(sum((V*L'*B).^2,2));
lam = (1/(max(tmp)*0.01));%^(1/p);                            %alpha
kesi = mean(sum(U.^2,2))*ones(size(V,1),1);%1e-8*ones(size(V,1),1);%
gamma = 1e17*ones(size(V,1),1);
rou_old = rou;
S_old = S;


% approximation of inverse
Lambda_MAX = eigs(V'*V,1,'lm');
LLt = L*L'; LtB = L'*B;
mu = 0.9/(rou*Lambda_MAX);
Precision = mu*speye(Nsource) - mu^2*L'/(eye(Nsensor) + mu*LLt)*L;

tic
w = 1;
U_inter_old = zeros(size(V,1),T); 
for iter = 1 : Max_iter 
%% ADMM (Source Update)
alpha = 0.6;
   for iter_ADMM = 1:ADMM_iter
%--------------------S update--------------------%
        S = Precision*(LtB + (1/mu)*S - rou*V'*(V*S - U + Z));
%--------------------U update--------------------%
        VS = V * S;
        VS_hat = alpha*VS + (1-alpha)*U_old;
        w = 1;
        U_inter_old = zeros(size(V,1),T); 
        for iter_U = 1:10%U_iter
% ============================= reweigthed L1 norm optimization ==============================%                
%            U = proxl21ARD(VS_hat + Z, w.^2, lam^(-p), rou);
%            rr = sqrt(kesi+sum(U.^2,2));
% %            w = rr.^(p-2).*sqrt(sum(U.^2,2));
% % %              w = p*rr.^(p-2).*sqrt(sum(U.^2,2));
% w = p*lam^(-p)*sqrt(sum(U.^2,2)+kesi).^(p-1);
%            w(w<max(w)*1e-6) = max(w)*1e-6;
% ============================= reweigthed L2 norm optimization ==============================%                      
            U = repmat(1./(rou/2 +2*w),1,size(U,2)).*(VS_hat + Z)*rou/2;
            rr = sqrt(kesi+sum(U.^2,2));
%             w = lam^(-p)*rr.^(p-2);
            w = p*lam^(-p)*rr.^(p-2);
            w(w<max(w)*1e-6) = max(w)*1e-6;
% ============================================================================================%            
            if  norm((U - U_inter_old),'fro')/norm(U,'fro') < 1e-6 || sum(w) == 0
                break;
            end   
            U_inter_old = U;
        end
%             rr = sqrt(kesi+sum(U.^2,2));
%             w = lam^(-p)*rr.^(p-2);
% %             w = p*lam^(-p)*rr.^(p-2);
%             U = repmat(1./(rou/2 +2*w),1,size(U,2)).*(VS_hat + Z)*rou/2;

%--------------------Z update--------------------%
        Z = Z + (VS_hat - U);
%--------------------stop criterion--------------------%
        primerror = norm(VS - U,'fro');
        dualerror = norm(rou*V'*(U - U_old),'fro');
        U_old = U;
        
        tol_prim = 1e-6*max(norm(U,'fro'),norm(VS,'fro'));
        tol_dual = 1e-6*rou*norm(V'*Z,'fro');
        
        Rprimerror = primerror/max(norm(U,'fro'),norm(VS,'fro'));
        Rdualerror = dualerror/norm(rou*V'*Z,'fro');
        
        if ~QUIET && mod(iter_ADMM,10) == 0
            fprintf('ADMM_wL21 : iteration: %g, Data-Fit : %g, PrimError: %g, DualError: %g\n', ier_ADMM, norm(B - L*S,'fro')/norm(B,'fro'), primerror, dualerror);
        end
        
        if primerror < tol_prim && dualerror < tol_dual
            break;
        end
        
%--------------------rou update--------------------%
        if rou_update && mod(iter_ADMM,10) == 0
           
            
            ratio = -1;
            if Rdualerror~=0
                ratio = sqrt(Rprimerror/Rdualerror);
            end
            
            tau_max = 2;
            if ratio>=1 && ratio<tau_max, tau = ratio;
            elseif ratio>1/tau_max && ratio<1, tau = 1/ratio;
            else tau = tau_max;
            end
            
            if Rprimerror > 10 * Rdualerror
                rou = tau*rou; Z = Z./tau;
            elseif Rdualerror > 10 * Rprimerror
                rou = rou/tau; Z = Z.*tau;
            end
            if ~QUIET
                fprintf('rou = %g, Rprimerror = %g, Rdualerror = %g\n',rou, Rprimerror, Rdualerror);
            end
            if rou ~= rou_old
                mu = 0.9/(rou*Lambda_MAX);
                Precision = mu*speye(Nsource) - mu^2*L'/(eye(Nsensor) + mu*LLt)*L;
            end
            rou_old = rou;
        end
    end

%% Gamma, Kesi and Alpha update
VS = V*S;
rr = sqrt(sum(VS.^2,2)+kesi);
gamma = p*lam^(-p)*rr.^(p-2);

[kesi,~,ldT] = diaginv_lanczos(L,1,V,gamma,250,1);
kesi = T * kesi;
fprintf('Maxgamma: %g,  ratio: %g, Maxkesi: %g, kesiratio: %g\n',max(1./gamma),max(1./gamma)/min(1./gamma),max(kesi),max(kesi)/min(kesi));%

rr = sqrt(sum(VS.^2,2)+kesi);


cost(iter+1) = size(B,2) * ldT  - sum((p-2)*lam^(-p)*rr.^p) + 2*size(V,1)*log(lam) + norm(B - L*S,'fro')^2 + trace((VS)'.*repmat(gamma',size(B,2),1)*(VS));
MSE = abs(cost(iter+1) - cost(iter))/abs(cost(iter+1));%norm(S - S_old,'fro')/norm(S,'fro');
fprintf('ADMM-VSSI-GGD : iteration =  %g, alpha =  %g, rou =  %g， MSE =  %g,\n',iter,lam,rou,MSE);
SS{iter} = S;
S_old = S;
% if abs(cost(iter+1) - cost(iter))/abs(cost(iter+1)) < tol, break; end
if MSE < tol, break; end



toc
end
% figure
% plot(cost(2:end));
S = SS;
end


function V = VariationEdge(VertConn)
Nsource = size(VertConn,1);
Nedge = numel(find(VertConn(:)~=0))/2;
V = sparse(Nedge,Nsource);
edge = 0;
for i = 1:Nsource
    idx = find(VertConn(i,:)~=0);
    idx = idx(idx<i);
    for j = 1:numel(idx)
        V(j+edge,i) = 1;
        V(j+edge,idx(j)) = -1;
    end
    edge = edge + numel(idx);
end
end

function Z = proxl21ARD(Y,w,lam,rou)
 [m,n] = size(Y);
 temp = lam*sqrt(w)./sqrt(sum(Y.^2,2))/rou;
 scale = pos(ones(m,1) - temp);
 Z = Y.*repmat(scale,1,n);
end

function x = pos(y)
x = max(0,y);
end
