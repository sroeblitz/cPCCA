function [chi,A,fopt]=opt_soft(EVS,flag, solver)
% computation of membership matrix by PCCA+
%
% [chi,A,fopt]=opt_soft(EVS,flag, solver)
%
% input parameter:
%    EVS:        (N,k)-matrix with eigenvectors (column by column)
%    flag:       0 = unscaled initial guess
%                1 = full optimization
%    solver:     solver for unconstrained optimization problem;
%                one of the following can be chosen:       
%                      'Nelder-Mead'
%                      'Levenberg-Marquardt'
%                      'Gauss-Newton'
% 
% output:
%    chi:        (N,k)-matrix with membership values
%    A:          (k,k)-tranformation matrix A s.th. chi=EVS*A
%    fopt        optimal function value  
%
% cite:
%
% [1] P. Deuflhard, M. Weber: Robust Perron Cluster Analysis in Conformation
%     Dynamics. Lin. Alg. App. 2005, 398c, Special issue on matrices and
%     mathematical biology, 161-184.
% [2] S. Roeblitz and M. Weber: Fuzzy Spectral Clustering by PCCA+. 
%     In: Classification and Clustering: Models, Software and Applications, 
%     WIAS Report No. 26, Berlin 2009.
%
% written by
% Susanna Roeblitz and Marcus Weber, Zuse Institute Berlin, Takustrasse 7, 14195 Berlin
% updated version: August, 2012


[N,k]=size(EVS);

if (k > 1)
        
     
    % search start simplex vertices ('inner simplex algorithm')
    index=indexsearch(EVS);
    
    
    % compute transformation matrix A as initial guess for local
    % optimization (maybe not feasible) 
    A=EVS(index,:);
    A=inv(A);
       
    % reduce optimization problem to size (k-1)^2
    alpha=zeros(1,(k-1)^2);
    for i=1:k-1
        for j=1:k-1
            alpha(j + (i-1)*(k-1)) = A(i+1,j+1);
        end
    end
    fopt=objective(alpha,EVS);
        
    if (flag > 0)   
        % perform optimization
        switch lower(solver)
            case 'nelder-mead'
                options = optimset('maxiter',2000,'TolFun',1e-8,'TolX',1e-8);
                [alpha,fopt] = fminsearch(@(alpha) objective(alpha,EVS),alpha,options);
                %disp(output.message)
            case 'levenberg-marquardt'    
                options = optimset('Algorithm','levenberg-marquardt','maxiter',500,'TolFun',1e-8,'TolX',1e-8);
                [alpha,~,fopt] = lsqnonlin(@(alpha) objective(alpha,EVS),alpha,[],[],options);
            case 'gauss-newton'
                if k>2
                    par.evs = EVS;
                    problem = 'problem_pcca_nlscon';
                    alpha = main_nlscon(alpha, par, problem);
                    fopt = feval(problem,alpha, '', par) ;
                else
                    disp('Gauss-Newton method does not work for k=2.')
                    disp('Instead, the Nelder-Mead algorithm is used.')
                    options = optimset('maxiter',2000,'TolFun',1e-8,'TolX',1e-8);
                    [alpha,fopt] = fminsearch(@(alpha) objective(alpha,EVS),alpha,options);
                end
                
            otherwise
                disp('Unknown solver for optimization problem.')
        end
        % complete A to meet constraints (positivity, partition of unity)
        for i=1:k-1
            for j=1:k-1
                A(i+1,j+1)=alpha(j + (i-1)*(k-1));
            end
        end
        A=fillA(A, EVS);
    %else
    %    fopt = -1;
    end
    

    chi=EVS*A;

else % special case: k==1
    %if (flag==1)
      fopt = 1.0;
    %else
    %  fopt = -1.0;
    %end
    chi=ones(N,1);
    
end

  
