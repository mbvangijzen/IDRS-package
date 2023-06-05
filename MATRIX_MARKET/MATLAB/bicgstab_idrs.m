%
% Comparison of IDRS-implementation of BiCGSTAB (Matlab and fortran) with
% Matlab's build-in implementation. The example is the (single frequency)
% wedge problem. Grids can be chosen for frequencies 1Hz to 6Hz.
%
% This driver needs the following files: idrs.m, P_mul.m, comp_om.m, make_P.m,
%                                        mmread.m (matrix-market routine)
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

clear all;
close all;
clc;

pathname = '../WEDGE_DATA/';

Title = 'Parameters';
example  = 'wedge';
s = '1';
grid = '1';
tol = '1e-6';
maxit = '10000';

contin = 1;
while contin
%% Menu
   prompt = { 'Grid (1-6):', ...
              'Tolerance:',...
              'Maximum iterations '};
   defaults = { grid, tol, maxit};
   lineNo=[1 50];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) contin = 0; break; end;
   grid     = char(params(1));
   tol      = char(params(2));
   maxit    = char(params(3));

%
%    Initialize figures
%  disp(' ');
   scrsz = get(0,'ScreenSize');
   fig1 = figure('Position',[scrsz(1) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);
   title('Convergence');
   xlabel('Number of iterations');
   ylabel('Residual norm');
   grid on;
   hold on;

%    Output: description of experiment:
   disp(['Problem = wedge', grid]);
%
%       Read the matrices:
   matrix = [pathname, example ,grid,'.mtx'];
   rhs = [pathname, example ,grid,'_b.mtx'];
   [A, n, m, entries, rep, field, symm] = mmread(matrix);
   [rhs, n, n_rhs, entries, rep, field, symm] = mmread(rhs);

%
   b = rhs(:,1);
%
   t0 = cputime;
   method = 2;
   [x,flag,relres,iter,resvec] = ...
      idrs(A,b,[],str2num(s),str2num(tol),str2num(maxit),2);
   t1 = cputime-t0;

   disp(' ');
   method1 = ['BiCGSTAB (Matlab)'];
   disp(['Method                  = ',method1]);
   disp(['Number of iteration     = ',num2str(iter)]);
   disp(['Relative residual norm  = ',num2str(relres)]);
   if ( flag > 0 )
      if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
      if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
      if ( flag == 3 ) disp('Break down!'); end;
   end
   disp(['CPU time = ', num2str(t1) ]);

% Plot convergence
   method1 = cellstr(method1);
   figure(fig1);
   x_as = [0:1:iter];
   plot(x_as,log10(resvec(1:iter+1))),'b';
%
   t0 = cputime;
   [x,flag,relres,iter,resvec] = ...
      bicgstab(A,b,str2num(tol),str2num(maxit));
   t1 = cputime-t0;

   disp(' ');
   method2 = ['BiCGSTAB (build-in)'];
   disp(['Method                  = ', method2]);
   disp(['Number of iteration     = ',num2str(iter)]);
   disp(['Relative residual norm  = ',num2str(relres)]);
   if ( flag > 0 )
      if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
      if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
      if ( flag == 3 ) disp('Stagnation'); end;
   end
   disp(['CPU time = ', num2str(t1) ]);
   method2 = cellstr(method2);

% Plot convergence
   figure(fig1);
   x_as = [0:1:size(resvec)-1];
   plot(x_as,log10(resvec),'r');

   disp(' ');
   method3 = ['BiCGSTAB (Fortran)'];
   matrix = [pathname, example ,grid,'.mtx'];
   run_mm_idrs = ['!../FORTRAN/mm_idrs -matrix ', matrix, ...
                     ' -bicgstab ', ...
                     ' -tol ', tol, ' -maxit ', maxit, ' -conv'];
   eval(run_mm_idrs);
   convergence = [pathname, example ,grid,'_c.mtx'];
   [resvec, iter, nprocs, entries, rep, field, symm] = mmread(convergence);
   method3 = cellstr(method3);

% Plot convergence
   figure(fig1);
   x_as = [0:1:size(resvec)-1];
   plot(x_as,log10(resvec),'k');
%
   legend([method1,method2,method3]);
end
