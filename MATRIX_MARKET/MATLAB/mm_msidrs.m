%
% Generic driver to test the (shifted) msidrs-algorithms on matrix-market matrices.
% To test the algorithms select a matrix-market matrix.
% Files with _b in the name correspond to the rhs-vector(s) and should not be selected.
% If no rhs-file (corresponding to the matrix-file) exists, % the vector b = 1 is taken as rhs.
% Files with _s in the name correspond to the shifts and should not be selected.
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

clear all;
close all;

clear all; close all;
[filename, pathname] = uigetfile('*.mtx','Choose a matrix');
[dum, matrix] = fileparts(filename);
file = [pathname, filename];

rhsname   = [pathname, matrix, '_b.mtx'];
shiftname = [pathname, matrix, '_s.mtx'];

tol = 1e-6;
max_it = 4000;

contin = 1;

variant = 'msidrs';
s = '4';
tol = '1e-6';
maxit = '10000';
in_s = '4';
in_tol = '1e-1';
in_it = '0';
language = 'Fortran';

while ( contin )
%% Menu
   Title = 'MSIDRS parameters';
   prompt = { 'Method: msidrs/msbicgstab/msqmridr ',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:'...
              'Omega shifts (y/n):'...
              'Language: (Matlab/Fortran)'};
   defaults = { variant, s, tol, maxit, in_s, in_tol, in_it, language };
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
   language = char(params(8));

   method = 0;
   if ( strcmpi(variant,'msidrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'msbicgstab' ) )
      method = 2;
      s = 1;
   end
   disp(['Language                = ',language]);

   if ( strcmpi(language,'fortran' ) )
      run_mm_msidrs = ['!../FORTRAN/mm_msidrs -matrix ', file, ...
                    ' -', variant, ' ', s, ...
                    ' -tol ', tol, ' -maxit ', maxit, ' -in_s ', in_s, ...
                    ' -in_tol ', in_tol, ' -in_it ', in_it];
      eval(run_mm_msidrs);
   else
%
% Read the matrix, shifts and rhs-vector:
      [A, n, m, entries, rep, field, symm] = mmread(file);
      [shift, n, n_shift, entries, rep, field, symm] = mmread(shiftname);
      if exist(rhsname) 
         [rhs, n, m, entries, rep, field, symm] = mmread(rhsname);
         b = rhs(:,1);
      else
         b = ones(n,1);
      end

%
% Output: description of experiment:
      disp(['Testproblem = ', matrix]);
      disp(['Solution method is ', variant, ' with s = ', s]);

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
      disp(' ');
      disp(['Number of iteration     = ',num2str(iter)]);
      for i = 1:n_shift
         disp(['Relative residual norm for shift ', num2str(shift(i)), ' = ',num2str(relres(i))]);
      end
      if ( flag > 0 ) 
         if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
         if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
         if ( flag == 3 ) disp('Break down!'); end;
      end
      disp(' ');
      disp(['CPU time                = ',num2str(time),'s.']);
      disp('==================================================');
   end
end
