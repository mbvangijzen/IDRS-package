%
% IDRS parameterised testproblems
%
% The software is distributed without any warranty.
%
% Martin van Gijzen
% Copyright (c) December 2023
%

clear all;
close all;
clc;

disp('FDM discretisation of a 3D convection-diffusion-reaction problem on a unit cube');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');

% Define system

% Defaults:
m = 64;
eps = 0.005;
beta(1) = 1;
beta(2) = 1;
beta(3) = 1;

% Defaults for the iterative solvers:

tol = 1e-6;
maxit = 1000;

Title = 'Parameters';
refine = 1;
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
   prompt = { 'Refinement in each direction:', ...
              'Method: idrs/bicgstab/qmridr',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:'};
   defaults = { num2str(refine), variant, num2str(s), num2str(tol), num2str(maxit), num2str(in_s), num2str(in_tol), num2str(in_it)};
   lineNo=[1 50];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) contin = 0; break; end;
   refine   = str2num(char(params(1)));
   variant  =         char(params(2));
   s        = str2num(char(params(3)));
   tol      = str2num(char(params(4)));
   maxit    = str2num(char(params(5)));
   in_s     = str2num(char(params(6)));
   in_tol   = str2num(char(params(7)));
   in_it    = str2num(char(params(8)));

   method = 0;
   if ( strcmpi(variant,'idrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'bicgstab' ) )
      method = 2;
      s = 1;
   end

% Generate matrix
   m = 64*refine;
   h = 1/(m+1);

   n = m*m*m;
   Sx = gallery('tridiag',m,-eps/h^2-beta(1)/(2*h),2*eps/h^2,-eps/h^2+beta(1)/(2*h));
   Sy = gallery('tridiag',m,-eps/h^2-beta(2)/(2*h),2*eps/h^2,-eps/h^2+beta(2)/(2*h));
   Sz = gallery('tridiag',m,-eps/h^2-beta(3)/(2*h),2*eps/h^2,-eps/h^2+beta(3)/(2*h));
   Is = speye(m,m);
   A = kron(kron(Is,Is),Sx) + kron(kron(Is,Sy),Is)+ kron(kron(Sz,Is),Is);

   x = linspace(h,1-h,m);
   sol = kron(kron(x.*(1-x),x.*(1-x)),x.*(1-x))';
   sol = sqrt(sol);
   b = A*sol;

   n_shift = 6;
   shift = [0, 1, 2, 3, 4, 5];

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
%
   t = cputime;
   if ( method == 1 || method == 2 )
      [x,flag,relres,iter,resvec] = ...
          msidrs(A,b,shift,s,[],tol,maxit,method);
   else
      [x,flag,relres,iter,resvec] = ...
          msqmridr(A,b,shift,s,[],tol,maxit,in_s,in_tol,in_it);
   end

   time = cputime - t;
   disp(['Multi-shift test problem, size is ',num2str(n)]);
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
