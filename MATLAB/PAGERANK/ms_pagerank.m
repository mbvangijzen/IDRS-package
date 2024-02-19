%
% This programme illustrates the use of the idrs-routines on an academic pagerank problem.
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

clear all;
close all;
clc;

disp('           Multishift pagerank problem');
disp('++++++++++++++++++++++++++++++++++++++++++++++++++');

maxlinks = 10;

% Defaults for the iterative solvers:

tol = 1e-6;
maxit = 1000;

Title = 'Parameters';
N = 1000;
p = [0.85, 0.95];
variant = 'idrs';
s = 4;
tol = 1e-6;
maxit = 10000;
in_s = 4;
in_tol = 1e-1;
in_it = 0;

contin = 1;
while contin
%% Menu
   prompt = { 'Nodes:', ...
              'Damping factor:',...
              'Method: idrs/bicgstab/qmridr',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:'};
   defaults = { num2str(N), num2str(p), variant, num2str(s), num2str(tol), num2str(maxit), num2str(in_s), num2str(in_tol), num2str(in_it)};
   lineNo=[1 50];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) contin = 0; break; end;
   N        = str2num(char(params(1)));
   p        = str2num(char(params(2)));
   variant  =         char(params(3));
   s        = str2num(char(params(4)));
   tol      = str2num(char(params(5)));
   maxit    = str2num(char(params(6)));
   in_s     = str2num(char(params(7)));
   in_tol   = str2num(char(params(8)));
   in_it    = str2num(char(params(9)));

   method = 0;
   if ( strcmpi(variant,'idrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'bicgstab' ) )
      method = 2;
      s = 1;
   end

% Generate matrix
   nnz_row = randi(maxlinks,N,1);
   nnz = sum( nnz_row );
   rows = zeros(nnz,1);
   cols = zeros(nnz,1);
   k = 0;
   for i = 1:N
      cols(k+1:k+nnz_row(i)) = randi(N,nnz_row(i),1);
      rows(k+1:k+nnz_row(i)) = i;
      k = k + nnz_row(i);
   end
   G = sparse(rows,cols,1,N,N);
   outlinks = sum(G,1)';
   ind = find( outlinks == 0 );
   outlinks(ind) = 1;
   D = spdiags(outlinks,0,N,N);
   A = G;
   shift = 1./p;
   n_shift = length(shift);
   b = ones(N,1);

%    Initialize figures
   scrsz = get(0,'ScreenSize');
   fig1 = figure('Position',[scrsz(1) 0.5*scrsz(4) 0.4*scrsz(3) 0.5*scrsz(4)]);
   if ( strcmpi(variant,'bicgstab') )
      title(['Convergence bicgstab, ']);
   else
      title(['Convergence ', variant, s]);
   end
   xlabel('Number of iterations');
   ylabel('Scaled residual norm');
   grid on;
   hold on;
%
   t = cputime;
   if ( method == 1 || method == 2 )
      [x,flag,relres,iter,resvec] = ...
          msidrs(A,b,shift,s,D,tol,maxit,method);
   else
      [x,flag,relres,iter,resvec] = ...
          msqmridr(A,b,shift,s,D,tol,maxit,in_s,in_tol,in_it);
   end
   x = D*x;
   for i_shift = 1:n_shift
      x(:,i_shift) = x(:,i_shift)/sum(x(:,i_shift));
   end

%
%    Output: description of experiment:
%
   time = cputime - t;
   disp(['Pagerank test problem, size is ',num2str(N)]);
   disp(['Solution method is ', variant, ' with s = ', num2str(s)]);
   disp(['Elapsed time            = ',num2str(time),'s.']);
   disp(['Number of iteration     = ',num2str(iter)]);
   for i = 1:n_shift 
      disp(['Relative residual norm for shift ', num2str(shift(i)), ' = ', num2str(relres(i))]);
   end
   if ( flag > 0 )
      if ( flag == 1 ) disp('Maximum number of iterations reached!'); end;
      if ( flag == 2 ) disp('Accuracy above prescribed tolerance!'); end;
      if ( flag == 3 ) disp('Break down!'); end;
   end
%
   disp('==================================================');

%
% Plot convergence
   figure(fig1);
   for i_shift = 1:n_shift
      x_as = [0:1:iter];
      plot(x_as,log10(resvec(:,i_shift)/resvec(1,i_shift)));
%     legend(['shift=',num2str(shift(i_shift))]);
   end

end
