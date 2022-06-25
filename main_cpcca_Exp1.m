%% Main run file Example 1 

close all; clear all;
tic

load('Exp1.mat')
n=size(P,1);

%density for normalization of eigenvectors
ovec=ones(n,1)/n; % (!)

%input parameters for PCCA+
target=1+eps;           % target eigenvalue: target=1.0000
kmin=3;                 % minimum number of clusters
kmax=3;                 % maximum number of clusters
solver2='gauss-newton'; % optimization method for optimal cluster number
opts.disp=0;            % no display of eigenvalue solver

% the following two parameters are not used if kmin==kmax
flag=1;                 % flag=1: full optimization, flag=0: unscaled initial guess
solver1='nelder-mead';  % optmization method for finding optimal cluster number

%Compute sub-space
[EVS,la,kmax]=compute_subspace(P,kmax,target,opts); 
la=diag(la);
%preprocess eigenvectors: separate real and imaginary parts
EVS=preprocessEVS(EVS,la);

% Optimization step by PCCA+: Compute A such that chi=EVS*A 
[chi,A,optval,EVS]=pcca(EVS,ovec,kmin,kmax,flag,solver1,solver2);

%% Sort chi vectors columnwise
[~,m]=max(chi);
[i,idx]=sort(m);
chi=chi(:,idx);

fprintf('\nList of states and their cluster with highest membership:\n')
for k=1:n
    [~,idx]=max(chi(k,:));
    fprintf('%d: %d\n',k,idx)
end


disp (' ')
disp ('Coarse-grainded transition matrix:') 
%Coarse-grained matrix:
Pc=inv(chi'*diag(ovec)*chi)*chi'*diag(ovec)*P*chi  

toc
et=toc;

% plot membership functions chi 
figure(1)
plot(chi,'-s','LineWidth',4,'Markersize',20)
set(gca,'FontSize',20)
%title('Membership functions chi','FontSize',20)
xlabel('States','Fontsize',20)
ylabel('Membership to cluster in range [0  1]','FontSize',20)
axis([1 6 0 1])
xticks([1 2 3 4 5 6])
xticklabels({'1','2','3','4','5','6'})
