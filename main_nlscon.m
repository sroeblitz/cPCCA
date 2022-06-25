function [final_par] = main_nlscon(xguess, par, problem)
%
% input to nlscon driver routine:
% -------------------------------
% 	xguess - initial parameter guess
% 	par - measurement timepoints and values
% 	problem  - file name of model description 
%                   (rhs evaluation and jacobian calculation)
%
% ouput:
% ------
%       final_par - resulting vector of iteration run
%
%     ____________________________________________________________ 
%
%
%
%*  Written by        L. Weimann 
%*  Purpose           Testexample for code NLSCON
%*  Version           2.3
%*  Revision          June 2004
%*  Latest Change     June 2004
%*  Library           CodeLib
%*  Code              Matlab 6.0
%*  Environment       Standard Matlab 6.0 environment on PC's,
%                     workstations and hosts.
%
%   modified by C. Stoetzel, T. Dierkes, 2011
%
%     ____________________________________________________________
%

  fidout = fopen(strcat('nlscon.m.out',num2str(xguess(1))),'w');
  fiddat = fopen(strcat('nlscon.m.dat',num2str(xguess(1))),'w');
  fprintf(1,'%s\n',' monitor: nlscon.m.out , data: nlscon.m.dat') ;
  
  xtol = 1.0e-8 ;
  
%     Execution mode: 0=Standard Mode, 1=Stepwise mode
  iopt.mode = 0 ;
  
%       Jacobian: 0=(same as value 3)
%                 1=supplied by user routine JAC
%                 2=computed by numerical differentation (no feedback) 
%                 3=computed by numerical differentation (with feedback)
  iopt.jacgen =3 ;
  
%     A posteriori statistical analysis: 0 = no, 1 = yes
  iopt.mprsta = 1 ;
  
%     Problem classification:
%     1 = linear , 2 = mildly nonlinear  3 = highly nonlinear
  iopt.nonlin = 3 ;
  
%     Broyden updates: 0 = inhibit, 1=allow
  iopt.qrank1 = 0 ;
  
%     Set output level on various units
  iopt.mprerr = 3 ;
  iopt.mprmon = 1 ;
  iopt.mprsol = 1 ;
  iopt.mprtim = 0 ;
  iopt.mprsta = 0 ;
  
%     Set output units
  iopt.luerr = fidout ;
  iopt.lumon = fidout ;
  iopt.lusol = fiddat ;
  iopt.lutim = fidout ;
  
%     Automatic row scaling of linear system (constraints)
%     and user scaling (measurements):
%     0 = allowed , 1 = inhibited
  iopt.norowscal = 0 ;

%

 wk.fcmin=1e-4;
%     ____________________________________________________________  
%
  
  x = xguess(:) ;
  fobs=0;
%   fobs = par(:,2) ;
%   tpoints = par(:,1);
  
    
%   fprintf(fidout,' Simulated  experimential  data  :%s\n',' ') ;
%   fprintf(fidout,' %10.7f\n',fobs) ;
%   fprintf(fidout,'%s\n\n\n',' ') ;
  xscal = ones(size(x)) ;
  fscal = ones(size(fobs)) ;
  
   

  f = feval(problem,x,'',par) ;
  m = length(f) ;
  
  info.ierr = -1 ;
  wk = [] ;
%   i = 0 ;
  stime = cputime ;
  
%   EVS=par.evs;
%   k=par.k;
%   while info.ierr == -1 
    [x,info,wk] = nlscon(m,problem,x,xscal,fobs,fscal,xtol,iopt,par,wk,info);
%     %%%%%%%%%%%%%%%%%%%%%
%     A=zeros(k,k);
%     % Matrix aus x rekonstruieren
%     for i=1:k-1
%         A(i+1,2:k)=x(((i-1)*(k-1)+1):i*(k-1));
%     end
%     [A,~]=fillA(A, EVS);
%     %chi=EVS*A;
%     F=trace(diag(1./A(1,:))*(A'*A));
%     f=k-F;
%     %%%%%%%%%%%%%%%%%%%%%
%     i = i+1;
% 	fprintf('Returned from call %i of NLSCON\n',i) ;
%   end
  etime = cputime ;
  cptime = etime-stime ;
  if cptime ~= 0
	fprintf(fidout, 'Time used = %9.3f Sec\n\n',cptime) ;
	fprintf('Time used = %9.3f Sec\n\n',cptime) ;
  end
 
  final_par = x(:) ;
   
  fclose(fidout);
  fclose(fiddat);
  
