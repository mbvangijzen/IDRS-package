%
% Multi-frequency Wedge test problems. 
%
% Files are provided on six grids for increasinly high frequencies:
% frequencies 1, 2, 3, 4, 5, and 6 Hz.
% Boundary conditions are Neumann (refecting) or Sommerfeld (absorbing)
%
% The data files are supplied in Matrix Market format
%
% The wedge problem is a classical test problem form geophysics.
% The formulation as a sequence of shifted problems can be found in:
%    M. Baumann and M.B. van Gijzen. 
%    Nested Krylov Methods for Shifted Linear Systems. 
%    SIAM Journal on Scientific Computing, 37(5):S90-S112, 2015. 
%
% This driver needs the following files: msidrs.m, msqmridr,
%                                        P_mul.m, comp_mu, make_P.m
%                                        mmread.m (matrix-market routine)
% First run setpath to be able to access these files.
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
grid = '4';
bc = 'neumann';
variant = 'msidrs';
s = '4';
tol = '1e-6';
maxit = '10000';
in_s = '4';
in_tol = '1e-1';
in_it = '0';

contin = 1;
while contin
%% Menu
   prompt = { 'Grid:', ...
              'Boundary condition: neumann/sommerfeld',...
              'Method: msidrs/msbicgstab/msqmridr',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:'};
   defaults = { grid, bc, variant, s, tol, maxit, in_s, in_tol, in_it};
   lineNo=[1 50];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) contin = 0; break; end;
   grid     = char(params(1));
   bc       = char(params(2));
   variant  = char(params(3));
   s        = char(params(4));
   tol      = char(params(5));
   maxit    = char(params(6));
   in_s     = char(params(7));
   in_tol   = char(params(8));
   in_it    = char(params(9));

   method = 0;
   if ( strcmpi(variant,'msidrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'msbicgstab' ) )
      method = 2;
      s = '1';
   end

   freq = zeros(str2num(grid),1);
   freq(1) = 1;
   for i = 2:str2num(grid)
      freq(i) = freq(i-1)*2;
   end

   if ( strcmpi(bc,'neumann') )
      pathname = '../../DATA/WEDGE_DATA/MULTI_NEUMANN/';
      shift = (2*pi*freq).^2;
   else
      pathname = '../../DATA/WEDGE_DATA/MULTI_SOMMERFELD/';
      shift = 2*pi*freq*sqrt(-1);
   end
   n_shift = length(shift);

%
%    Output: description of experiment:
%
%    Initialize figures
   scrsz = get(0,'ScreenSize');
   fig1 = figure('Position',[scrsz(1) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);
   if ( strcmpi(variant,'bicgstab') )
      title(['Convergence msbicgstab, ']);
   else
      title(['Convergence ', variant, s]);
   end
   xlabel('Number of iterations');
   ylabel('Scaled residual norm');
   grid on;
   hold on;
   fig2 = figure('Position',[0.6*scrsz(3) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);

%       Read the matrices:
   K_matrix = [pathname, 'wedge',grid,'_K.mtx'];
   M_matrix = [pathname, 'wedge',grid,'_M.mtx'];
   rhs      = [pathname, 'wedge',grid,'_b.mtx'];
   [K, n, m, entries, rep, field, symm] = mmread(K_matrix);
   [M, n, m, entries, rep, field, symm] = mmread(M_matrix);
   [b, n, nrhs, entries, rep, field, symm] = mmread(rhs);
%
   t = cputime;
   if ( method == 1 || method == 2 )
      [x,flag,relres,iter,resvec] = ...
          msidrs(K,b,shift,str2num(s),M,str2num(tol),str2num(maxit),method);
   else
      [x,flag,relres,iter,resvec] = ...
          msqmridr(K,b,shift,str2num(s),M,str2num(tol),str2num(maxit),...
                   str2num(in_s),str2num(in_tol),str2num(in_it));
   end

   time = cputime - t;
   disp(['Multi-frequency test problem, size is ',num2str(n)]);
   disp(['Solution method is ', variant, ' with s = ', s]);
   disp(['Elapsed time            = ',num2str(time),'s.']);
   disp(['Number of iteration     = ',num2str(iter)]);
   for i = 1:n_shift
      disp(['Results for frequency  = ', num2str(freq(i))]);
      disp(['Relative residual norm = ', num2str(relres(i))]);
%     disp(['Solution at source     = ', num2str(x(3*2^(str2num(grid)-1)+1,i))]);
   end
   if ( flag > 0 )
      if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
      if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
      if ( flag == 3 ) disp('Break down!'); end;
   end
%
   disp('==================================================');

%
% Plot solution:
   for i_shift = 1:n_shift
      nx = 3*2^str2num(grid)+1;
      ny = 5*2^str2num(grid)+1;
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
%     legend(['shift=',num2str(shift(i_shift))]);
   end

end
