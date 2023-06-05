%
% Multi-frequency Wedge test problems. 
%
% Files are provided on six grids for increasinly high frequencies:
% frequencies 1, 2, 3, 4, 5, and 6 Hz.
% Boundary conditions are Neumann (refecting) or Sommerfeld (absorbing)
%
% Files have the following names and meanings:
%
%     wedge?_neumann.mtx:      the system matrix for the neumann problem
%     wedge?_neumann_s.mtx:    the shifts for the neumann problem
%     wedge?_neumann_b.mtx:    the rhs-vector for the neumann problem
%     wedge?_sommerfeld.mtx:   the system matrix for the neumann problem
%     wedge?_sommerfeld_s.mtx: the shifts for the sommerfeld problem
%     wedge?_sommerfeld_b.mtx: the rhs-vector for the sommerfeld problem
%
% The data files are supplied in Matrix Market format
%
% The wedge problem is a classical test problem form geophysics.
% The formulation as a sequence of shifted problems can be found in:
%    M. Baumann and M.B. van Gijzen. 
%    Nested Krylov Methods for Shifted Linear Systems. 
%    SIAM Journal on Scientific Computing, 37(5):S90-S112, 2015. 
%
% This driver needs the following files: idrs.m, P_mul.m, comp_mu, make_P.m
%                                        mmread.m (matrix-market routine)
%
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

clear all;
close all;
clc;

Title = 'Parameters';
freq = '4';
bc = 'sommerfeld';
variant = 'msidrs';
s = '4';
tol = '1e-6';
maxit = '10000';
in_s = '4';
in_tol = '1e-1';
in_it = '0';
language = 'Fortran';
pathname = '../WEDGE_DATA/';

contin = 1;
while contin
%% Menu
   prompt = { 'Max. frequency:', ...
              'Boundary condition: neumann/sommerfeld',...
              'Method: msidrs/msbicgstab/msqmridr',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:',...
              'Language: (Matlab/Fortran)'};
   defaults = { freq, bc, variant, s, tol, maxit, in_s, in_tol, in_it, language};
   lineNo=[1 50];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) contin = 0; break; end;
   freq     = char(params(1));
   bc       = char(params(2));
   variant  = char(params(3));
   s        = char(params(4));
   tol      = char(params(5));
   maxit    = char(params(6));
   in_s     = char(params(7));
   in_tol   = char(params(8));
   in_it    = char(params(9));
   language = char(params(10));

   method = 0;
   if ( strcmpi(variant,'msidrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'msbicgstab' ) )
      method = 2;
      s = '1';
   end

%
%    Output: description of experiment:
   disp(['Language                = ',language]);
%
%    Initialize figures
   scrsz = get(0,'ScreenSize');
   fig1 = figure('Position',[scrsz(1) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);
   if ( strcmpi(variant,'bicgstab') )
      title(['Convergence msbicgstab, ', language]);
   else
      title(['Convergence ', variant, s, ', ', language]);
   end
   xlabel('Number of iterations');
   ylabel('Scaled residual norm');
   grid on;
   hold on;
   fig2 = figure('Position',[0.6*scrsz(3) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);

   if ( strcmpi(language,'Fortran' ) )
      matrix = [pathname, 'wedge',freq,'_',bc,'.mtx'];
      run_mm_msidrs = ['!../FORTRAN/mm_msidrs -matrix ', matrix, ...
                     ' -', variant, ' ', s, ...
                     ' -tol ', tol, ' -maxit ', maxit, ' -in_s ', in_s, ...
                     ' -in_tol ', in_tol, ' -in_it ', in_it, ' -save -conv'];
      eval(run_mm_msidrs);
      solution = [pathname, 'wedge',freq,'_',bc,'_x.mtx'];
      convergence = [pathname, 'wedge',freq,'_',bc,'_c.mtx'];
      [x, n, n_shift, entries, rep, field, symm] = mmread(solution);
      [resvec, iter, dum, entries, rep, field, symm] = mmread(convergence);
      iter = iter-1;
   else
%       Read the matrices:
      matrix = [pathname, 'wedge',freq,'_',bc,'.mtx'];
      rhs    = [pathname, 'wedge',freq,'_',bc,'_b.mtx'];
      shifts = [pathname, 'wedge',freq,'_',bc,'_s.mtx'];
      [A, n, m, entries, rep, field, symm] = mmread(matrix);
      [b, n, nrhs, entries, rep, field, symm] = mmread(rhs);
      [shift, n1, n_shift, entries, rep, field, symm] = mmread(shifts);
%
      t = cputime;
      if ( method == 1 || method == 2 )
         [x,flag,relres,iter,resvec] = ...
             msidrs(A,b,shift,str2num(s),str2num(tol),str2num(maxit),method);
      else
         [x,flag,relres,iter,resvec] = ...
             msqmridr(A,b,shift,str2num(s),str2num(tol),str2num(maxit),...
                      str2num(in_s),str2num(in_tol),str2num(in_it));
      end

      time = cputime - t;
      disp(['Test problem = ', matrix]);
      disp(['Results for ', variant, ' with s = ', s]);
      disp(['Elapsed time            = ',num2str(time),'s.']);
      disp(['Number of iteration     = ',num2str(iter)]);
      for i = 1:length(shift)
         disp(['Relative residual norm for shift ', num2str(shift(i)), ' = ', num2str(relres(i))]);
      end
      if ( flag > 0 )
         if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
         if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
         if ( flag == 3 ) disp('Break down!'); end;
      end
%
      disp('==================================================');
   end

%
% Plot solution:
   for i_shift = 1:n_shift
      nx = 3*2^str2num(freq)+1;
      ny = 5*2^str2num(freq)+1;
      if strcmp( bc, 'sommerfeld' )
         y = reshape(x(n/2+1:n,i_shift),nx,ny);
      else
         y = reshape(x(:,i_shift),nx,ny);
      end
      figure(fig2);
      imagesc(real(y));
      drawnow;
      pause(0.5);
   end
%
% Plot convergence
   figure(fig1);
   for i_shift = 1:n_shift
      x_as = [0:1:iter];
      plot(x_as,log10(resvec(:,i_shift)/resvec(1,i_shift)));
%      legend(['shift=',str2num(shift(i_shift)]);
   end

end
