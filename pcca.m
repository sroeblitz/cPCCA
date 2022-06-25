function [chiopt,A,optval,EVS_orth]=pcca(EVS,pi,kmin,kmax,flag,solver1,solver2)
% This algorithm generates a fuzzy clustering such that the resulting
% membership functions are as crisp (characteristic) as possible.
%
%
% Input:
%    EVS:        (N,kmax+x)-matrix with columns as invariant subpsace
%    pi:         vector for orthogonalization (e.g. stationary density)
%    kmin:       minimum number of clusters
%    kmax:       maximum number of clusters
%    flag:       0 = unscaled initial guess
%                1 = full optimization
%    solver1:    optimization method used to find the optimal cluster number;
%                (only relevant if flag=1);
%                one of the following can be chosen:       
%                      'Nelder-Mead'
%                      'Levenberg-Marquardt'
%                      'Gauss-Newton'
%    solver2:	 optimization method used to find A for the optimal cluster number  
%                one of the following can be chosen:       
%                      'Nelder-Mead'
%                      'Levenberg-Marquardt'
%                      'Gauss-Newton'
% 
% Output:
%   chiopt:     (N,kopt)-membership matrix
%   A:          (k,k)-transformation matrix with chi=EVS*A (EVS eigenvectors of A)    
%   optval:     value of the objective function: val=k_opt-trace(S) -> min
%   EVS_orth:   orthogonalized eigenvectors
%
% Cite:
%
% [1] P.Deuflhard, M.Weber: Robust Perron Cluster Analysis in Conformation Dynamics,
%     Lin. Alg. App. 2005, 398c, 161-184.
% [2] S. Roeblitz and M. Weber: Fuzzy spectral clustering by PCCA+: Application to
%     Markov state models and data classification. Advances in Data Analysis and
%     Classification 7(2):147–179, 2013. doi: 10.1007/s11634-013-0134-6.

% Written by Susanna Roeblitz and Marcus Weber, Zuse Institute Berlin, 2012
% updated version: June, 2022
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% orthogonalize eigenvectors w.r.t. pi; make 1st column constant 1
% such that EVS'*diag(pi)*EVS=Identity
EVS_orth=orthogon(EVS, pi);    
disp (['Subspace error after orthogonalization: ' num2str(subspace(EVS_orth,EVS))])

% decide for "best" number of clusters (other criteria are possible):
% the optimal value for trace(S) can be at most kt
% therefore crispness=trace(S)/kt<1, but should be as large as possible

disp (' ')
disp ('Optimize crispness vs. number of clusters')
disp ('=============================================')

opt_vec=zeros(1,kmax-kmin+1);

% compute optimal solution for every desired number of clusters
if kmax>kmin
    for kt=kmin:kmax
        [~,~,val]=opt_soft(EVS_orth(:,1:kt), flag, solver1);    %call to optimization routine
        disp (['For ' int2str(kt) ' clusters : crispness = ' num2str((kt-val)/kt)])
        disp (' ')
        opt_vec(kt-kmin+1)=(kt-val)/kt;
    end
    
    
    figure(1)
    h=plot(kmin:kmax,opt_vec,'-x');
    title('Number cluster vs. crispness (Value should be as large as possible)');
    xlabel('number of clusters')
    ylabel('crispness')
    
    disp(' ')
    kopt = input('Your choice for the number of clusters: ');
    disp (' ')
    if isempty(kopt) || kopt>kmax || kopt<kmin
        disp('invalid number of clusters')
        [~,kopt] = max(opt_vec);
        kopt = kmin + kopt -1;
        disp (' ')
        disp (['Optimum found for ' int2str(kopt) ' clusters : crispness = ' num2str(opt_vec(kopt-kmin+1))]) 
        disp (' ')
    else
        disp (['User optimum: ' int2str(kopt) ' clusters : crispness = ' num2str(opt_vec(kopt-kmin+1))]) 
        disp (' ')
    end
else
    kopt=kmax;
end

% repeat call to opt_soft for optimal cluster number kopt
%[chiopt,A,optval]=opt_soft(EVS_orth(:,1:kopt), flag, solver);
[chiopt,A,optval]=opt_soft(EVS_orth(:,1:kopt), 1, solver2);

disp ('Statistical weights of these clusters') 
weights=chiopt'*pi;
for i=1:length(weights)
    disp (['For ' int2str(i) '-th cluster : Weight = ' num2str(weights(i))])
end


