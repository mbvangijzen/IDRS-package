%
% Generic driver to test the (unshifted) idrs-algorithms on matrix-market matrices.
% To test the algorithms select a matrix-market matrix. 
% Files with _b in the name correspond to the rhs-vector(s) and should not be selected. 
% If no rhs-file (corresponding to the matrix-file) exists, 
% the vector b = 1 is taken as rhs.
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

clear all; close all;
[filename, pathname] = uigetfile('*.mtx','Choose a matrix');
[dum, matrix] = fileparts(filename);
file = [pathname, filename];

rhsname   = [pathname, matrix, '_b.mtx'];
H_name    = [pathname, matrix, '_h.mtx'];

tol = 1e-6;
max_it = 4000;

contin = 1;

variant = 'idrs';
s = '4';
tol = '1e-6';
maxit = '1000';
in_s = '4';
in_tol = '1e-1';
in_it = '0';
nritz = '0';
precon = 'n';
recycle = 'n';
language = 'Fortran';

while ( contin )
%% Menu
   Title = 'IDRS parameters';
   prompt = { 'Method: idrs/bicgstab/qmridr',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:',...
              'Number of Ritzvalues:',...
              'Preconditioner (y/n):',...
              'Recycle (y/n):',...
              'Language: (Matlab/Fortran)'};
   defaults = { variant, s, tol, maxit, in_s, in_tol, in_it, nritz, precon, recycle, language};
   lineNo=[1 50];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) contin = 0; break; end;
   variant  = char(params(1));
   s        = char(params(2));
   tol      = char(params(3));
   maxit    = char(params(4));
   in_s     = char(params(5));
   in_tol   = char(params(6));
   in_it    = char(params(7));
   nritz    = char(params(8));
   precon   = char(params(9));
   recycle  = char(params(10));
   language = char(params(11));

   method = 0;
   if ( strcmpi(variant,'idrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'bicgstab' ) )
      method = 2;
      s = '1';
   end
   disp(['Language                = ',language]);

   if ( strcmpi(language,'fortran' ) )
      run_mm_idrs = ['!../FORTRAN/mm_idrs -matrix ', file, ...
                    ' -', variant, ' ', s, ...
                    ' -tol ', tol, ' -maxit ', maxit, ' -in_s ', in_s, ...
                    ' -hessenberg ', nritz, ...
                    ' -in_tol ', in_tol, ' -in_it ', in_it, ' -save'];
      if ( strcmpi(precon,'y' ) )
         run_mm_idrs = [run_mm_idrs, ' -precon'];
      end
      if ( strcmpi(recycle,'y' ) )
         run_mm_idrs = [run_mm_idrs, ' -recycle'];
      end
      eval(run_mm_idrs);
      if ( str2num(nritz) > 0 & method == 1 )
         [H,n1, m1, entries, rep, field, symm] = mmread(H_name);
      end
   else
%
% Read the matrix and rhs-vector:
      [A, n, m, entries, rep, field, symm] = mmread(file);
      if exist(rhsname) 
         [rhs, n, n_rhs, entries, rep, field, symm] = mmread(rhsname);
      else
         rhs = ones(n,1);
      end
      if ( strcmpi(precon,'y' ) )
         D = diag(A);
         M = spdiags(D,0,n,n);
      else
         M = [];
      end

%
% Output: description of experiment:
      disp(['Testproblem = ', matrix]);
      disp(['Solution method is ', variant, ' with s = ', s]);

%
% Loop over different rhs-vectors
      x = zeros( n,n_rhs);
      x0 = zeros(n,1);
      t_tot = 0;
      U0 = zeros(n,str2num(s));
      for irhs = 1:n_rhs
%
% Copy rhs-vector:
         b = rhs(:,irhs);
%
         t0 = cputime;
         if ( method == 1 || method == 2 )
            if ( irhs > str2num(s) & strcmpi(recycle,'y' ) )
               U0 = x(:,irhs-str2num(s):irhs-1);
               [x(:,irhs),flag,relres,iter,resvec,H] = ...
                  idrs(A,b,M,str2num(s),str2num(tol),str2num(maxit),method,x0,U0);
            else
               [x(:,irhs),flag,relres,iter,resvec] = ...
                  idrs(A,b,M,str2num(s),str2num(tol),str2num(maxit),method,x0);
            end
            x0 = x(:,irhs);
         else
            [x(:,irhs),flag,relres,iter,resvec] = ...
                qmridr(A,b,M,str2num(s),str2num(tol),str2num(maxit),...
                    str2num(in_s),str2num(in_tol),str2num(in_it),x0);
                x0 = x(:,irhs);
         end
         t_tot = t_tot + cputime - t0;

         disp(' ');
         disp(['Results for system ',num2str(irhs)] );
         disp(['Number of iteration     = ',num2str(iter)]);
         disp(['Relative residual norm  = ',num2str(relres)]);
         if ( flag > 0 )
            if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
            if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
            if ( flag == 3 ) disp('Break down!'); end;
         end

      end
      disp(' ');
      disp(['CPU time                = ',num2str(t_tot),'s.']);
      disp('==================================================');
   end

   if ( str2num(nritz) > 0 & method == 1 )
      H = H(1:str2num(nritz),1:str2num(nritz));
      figure;
      plot(eig(H),'*');
      title('Ritzvalues');
      xlabel('Re');
      ylabel('Im');
   end
end
