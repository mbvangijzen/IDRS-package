%
% Ocean test problem. 
%
% The test problems are discretizations of Stommel's and Sag's ocean models. 
% Files are provided on six grids (spacing 1, 2, 3, 4, 5 and 6 degrees).
% The Stommel model on a 0.5 degree grid is also provided as stommel0.mtx
%
% Files have the following names and meanings:
%
%     stommel?.mtx:     the system matrix 
%     stommel?_b.mtx:   the right-hand sides. Every file contains 12 
%                       rhs vectors, one for each month.
%     stommel?_map.mtx: map from system dofs to earth dofs, needed to
%                       plot the solution
%     sag?.mtx:         the system matrix 
%     sag?_b.mtx:       the right-hand sides. Every file contains 12 
%                       rhs vectors, one for each month.
%     sag?_map.mtx: map from system dofs to earth dofs, needed to
%                       plot the solution
%     bathymetry.mtx:   contains waterdepth in each gridpoint (spacing 1 dg)
%
% The data files are supplied in Matrix Market format
%
% The discretisations are described in: 
%     M. B. van Gijzen, C. B. Vreugdenhil, and H. Oksuzoglu, 
%     The Finite Element Discretization for Stream-Function Problems 
%     on Multiply Connected Domains, 
%     J. Comp. Phys., 140, 1998, pp. 30-46. (copyright Academic Press)
% Please put a reference to this paper if you use thes test problems.
%
% This driver needs the following files: idrs.m, P_mul.m, comp_mu.m, make_P.m
%                                        mmread.m (matrix-market routine)
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

clear all;
close all;
clc;

month = char('January ','February','March', 'April', 'May', 'June', 'July', ...
             'August','September','October','November','December');

Title = 'Parameters';
example  = 'stommel';
grid = '4';
variant = 'idrs';
s = '4';
tol = '1e-6';
maxit = '10000';
in_s = '4';
in_tol = '1e-1';
in_it = '0';
language = 'Fortran';
pathname = '../OCEAN_DATA/';
[bathymetry, n1, m1, entries, rep, field, symm] = mmread([pathname,'bathymetry.mtx']);

contin = 1;
while contin
%% Menu
   prompt = { 'Example: stommel/sag',...
              'Grid (1-6):', ...
              'Method: idrs/bicgstab/qmridr',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:',...
              'Language: (Matlab/Fortran)'};
   defaults = { example, grid, variant, s, tol, maxit, in_s, in_tol, in_it, language};
   lineNo=[1 50];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) contin = 0; break; end;
   example  = char(params(1));
   grid     = char(params(2));
   variant  = char(params(3));
   s        = char(params(4));
   tol      = char(params(5));
   maxit    = char(params(6));
   in_s     = char(params(7));
   in_tol   = char(params(8));
   in_it    = char(params(9));
   language = char(params(10));

   method = 0;
   if ( strcmpi(variant,'idrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'bicgstab' ) )
      method = 2;
      s = '1';
   end
   disp(['Language                = ',language]);

%
%    Initialize figures
   scrsz = get(0,'ScreenSize');
   fig1 = figure('Position',[scrsz(1) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);
   if ( strcmpi(variant,'bicgstab') )
      title(['Convergence BiCGSTAB, ', language]);
   else
      title(['Convergence ', variant, s, ', ', language]);
   end
   xlabel('Number of iterations');
   ylabel('Residual norm');
   grid on;
   hold on;
   fig2 = figure('Position',[0.6*scrsz(3) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);

   if ( strcmpi(language,'Fortran' ) )
      matrix = [pathname, example ,grid,'.mtx'];
      run_mm_idrs = ['!../FORTRAN/mm_idrs -matrix ', matrix, ...
                     ' -', variant, ' ', s, ...
                     ' -tol ', tol, ' -maxit ', maxit, ' -in_s ', in_s, ...
                     ' -in_tol ', in_tol, ' -in_it ', in_it, ' -precon -save -conv'];
      eval(run_mm_idrs);
      solution = [pathname, example ,grid,'_x.mtx'];
      convergence = [pathname, example ,grid,'_c.mtx'];
      map = [pathname, example ,grid,'_map.mtx'];
      [x, n, n_rhs, entries, rep, field, symm] = mmread(solution);
      [resvec, iter, dum, entries, rep, field, symm] = mmread(convergence);
      iter = iter-1;
      [P, n1, m1, entries, rep, field, symm] = mmread(map);
      disp(' ');

% Plot convergence
      figure(fig1);
      x_as = [0:1:iter];
      plot(x_as,log10(resvec));
      legend(month);
      drawnow;

   else
%
% Read the matrices:
      matrix = [pathname, example, grid,'.mtx'];
      rhs = [pathname, example ,grid,'_b.mtx'];
      map = [pathname, example ,grid,'_map.mtx'];
      [A, n, m, entries, rep, field, symm] = mmread(matrix);
      [rhs, n, n_rhs, entries, rep, field, symm] = mmread(rhs);
      [P, n1, m1, entries, rep, field, symm] = mmread(map);
      D = diag(A);
      M = spdiags(D,0,n,n);

%
% Output: description of experiment:
      disp(['Testproblem = ', matrix]);
      disp(['Solution method is ', variant, ' with s = ', s]);
%
% Loop over different months
      x = zeros( n,n_rhs);
      x0 = zeros(n,1);
      t_tot = 0;
      for irhs = 1:n_rhs
%
% Copy rhs-vector:
         b = rhs(:,irhs);
%
         t0 = cputime;
         if ( method == 1 || method == 2 )
            [x(:,irhs),flag,relres,iter,resvec] = ...
                  idrs(A,b,M,str2num(s),str2num(tol),str2num(maxit),method,x0);
            x0 = x(:,irhs);
         else
            [x(:,irhs),flag,relres,iter,resvec] = ...
                qmridr(A,b,M,str2num(s),str2num(tol),str2num(maxit),...
                    str2num(in_s),str2num(in_tol),str2num(in_it),x0);
                x0 = x(:,irhs);
         end
         t_tot = t_tot + cputime - t0;

         disp(' ');
         disp(['Month is ',month(irhs,:)] );
         disp(['Number of iteration     = ',num2str(iter)]);
         disp(['Relative residual norm  = ',num2str(relres)]);
         if ( flag > 0 )
            if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
            if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
            if ( flag == 3 ) disp('Break down!'); end;
         end

% Plot convergence
         figure(fig1);
         x_as = [0:1:iter];
         plot(x_as,log10(resvec));
         legend(month(1:irhs,:));
         drawnow;
%
      end
      disp(' ');
      disp(['CPU time                = ',num2str(t_tot),'s.']);
      disp('==================================================');
   end

% Map solution on earth coordinates
   deg = str2num(grid);
   Px = P*x;
   if ( strcmp(example,'sag')  )
      Px = Px(1:n1/2,:);
   end
%
% Plot solution:
   for irhs = 1:n_rhs
      y = reshape(Px(:,irhs),360/deg+1,180/deg);
      figure(fig2);
      contour(0:deg:360,-(90-deg/2):deg:(90-deg/2),y',20)
      hold on;
      caxis([-1e5 5e4]);
      contour(0:360,-89.5:89.5,bathymetry',[0 0],'k')
      title(month(irhs,:));
      drawnow;
      pause(0.5);
      hold off;
   end
%

end
