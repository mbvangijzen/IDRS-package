%
% IDRS parameterised testproblems
%
% The software is distributed without any warranty.
%
% Martin van Gijzen
% Copyright (c) February 2014
%

clear all;
close all;

disp('FDM discretisation of a 3D convection-diffusion-reaction problem on a unit cube');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');

% Define system

% Defaults:
m = 64;
eps = 0.005;
beta(1) = 1;
beta(2) = 1;
beta(3) = 1;
r = 5;

Title = 'Problem Parameters';
prompt = {'Nodes in each direction:','Diffusion:','Convection x-direction:','Convection y-direction:','Convection z-direction:', 'Reaction:'};
defaults = {num2str(m),num2str(eps),num2str(beta(1)),num2str(beta(2)),num2str(beta(3)),num2str(r)};
lineNo=[1 25];
params=inputdlg(prompt,Title,lineNo,defaults);

% Actual values:
m = str2num(char(params(1)));
eps = str2num(char(params(2)));
beta(1) = str2num(char(params(3)));
beta(2) = str2num(char(params(4)));
beta(3) = str2num(char(params(5)));
r = str2num(char(params(6)));

% Generate matrix
h = 1/(m+1);

n = m*m*m;
Sx = gallery('tridiag',m,-eps/h^2-beta(1)/(2*h),2*eps/h^2,-eps/h^2+beta(1)/(2*h));
Sy = gallery('tridiag',m,-eps/h^2-beta(2)/(2*h),2*eps/h^2,-eps/h^2+beta(2)/(2*h));
Sz = gallery('tridiag',m,-eps/h^2-beta(3)/(2*h),2*eps/h^2,-eps/h^2+beta(3)/(2*h));
Is = speye(m,m);
I = speye(n,n);
A = kron(kron(Is,Is),Sx) + kron(kron(Is,Sy),Is)+ kron(kron(Sz,Is),Is) -r*I;

x = linspace(h,1-h,m);
sol = kron(kron(x.*(1-x),x.*(1-x)),x.*(1-x))';
sol = sqrt(sol);
b = A*sol;

disp(' ');
disp('The parameters of the problem are :');
disp(['Gridsize h = ',num2str(h),';']);
disp(['Number of equations = ', num2str(n),';']);
disp(['Diffusion parameter = ', num2str(eps),';']);
disp(['Convection parameters = (',num2str(beta(1)),',',num2str(beta(2)),',',num2str(beta(3)),');']);
disp(['Reaction parameter = ',num2str(r),';']);
disp(' ');

% Defaults for the iterative solvers:

fig = figure;
hold on;
xlabel('Number of MATVECS')
ylabel('log(|r|/|b|)')
title(['Instationary linear, \epsilon =', num2str(eps),'; \beta = ', num2str(beta),'; r =',num2str(r)]);
grid on;
c = ['b','g','r','c','m','y','k'];
kl = 0;
methods = [];

tol = 1e-6;
maxit = 1000;

method = 1;
s = 4;
dt = 1;
t_end = 10;
contin = 1;
while contin
   Title = 'Solver Parameters';
   prompt = {'IDRS/BiCGSTAB/QMRIDR (1/2/3):','s:','Initial guess (y/n):','Initial search space (y/n):',...
             'Ritz-omega (y/n):', 'End time:','Time step:'};
   defaults = {num2str(method),int2str(s),'y','n','n',num2str(t_end),num2str(dt)};
   lineNo=[1 25];
   params=inputdlg(prompt,Title,lineNo,defaults);
   if ( isempty(params) ) break; end;

% Actual values:
   method = str2num(char(params(1)));
   s = str2num(char(params(2)));
   if ( method == 2 ) s = 1; end;
   iniguess = ( char(params(3)) == 'y' );
   inispace = ( char(params(4)) == 'y' );
   user_omega = ( char(params(5)) == 'y' );
   t_end = str2num(char(params(6)));
   dt = str2num(char(params(7)));

% Initial value
   x = zeros(n,1);

% Euler backward time integration:
   t0 = cputime;
   J = speye(n,n) + dt*A;
   U0 = [];
   x0 = [];
   omega = [];
   t = dt;
   [x, flag, relres, iter, resvec,U0, omega] = idrs_ritz( J, dt*b, s, [], tol, maxit );
   disp(['Time: ', num2str(t), '; IDR(',num2str(s),') iterations: ',int2str(iter)]);
   disp(['Difference with stationary solution: ', num2str(norm(x-sol))]);
   tot_res = [resvec];

   if ( ~inispace ) U0 = []; end;
   if ( ~user_omega ) omega = []; end;

   while t < t_end
%
      if ( iniguess ) x0 = x; end;
      t = t + dt;
      if ( method == 1 ) 
         [x, flag, relres, iter, resvec] = idrs( J, (x+dt*b), s, [], tol, maxit, method, x0, U0, omega );
         disp(['Time: ', num2str(t), '; IDR(',int2str(s),') iterations: ',int2str(iter)]);
      elseif ( method == 2 ) 
         [x, flag, relres, iter, resvec] = idrs( J, (x+dt*b), 1, [], tol, maxit, method, x0, U0, omega );
         disp(['Time: ', num2str(t), '; BiCGSTAB iterations: ',int2str(iter)]);
      else
         [x, flag, relres, iter, resvec] = qmridr( J, (x+dt*b), s, [], tol, maxit, [], [], [], x0, omega );
         disp(['Time: ', num2str(t), '; QMRIDR(',int2str(s),') iterations: ',int2str(iter)]);
      end
      tot_res = [tot_res; resvec];
      disp(['Difference with stationary solution: ', num2str(norm(x-sol))]);
   end
   te = cputime;
   disp(['Total CPU time: ', num2str(te-t0)]);
   disp(' ');

   kl = kl + 1;
   if ( method == 1 )
      if ( inispace | user_omega )
         solver = ['Ritz-IDR(',int2str(s),') with recycling'];
      elseif ( user_omega )
         solver = ['Ritz-IDR(',int2str(s),')'];
      elseif ( inispace )
         solver = ['IDR(',int2str(s),') with recycling'];
      else
         solver = ['IDR(',int2str(s),')'];
      end
   elseif ( method == 2 )
      if ( inispace | user_omega )
         solver = 'Ritz-BiCGSTAB with recycling';
      elseif ( user_omega )
         solver = 'Ritz-BiCGSTAB';
      elseif ( inispace )
         solver = 'BiCGSTAB with recycling';
      else
         solver = 'BiCGSTAB';
      end
   else
      if ( user_omega )
         solver = ['Ritz-QMRIDR(',int2str(s),')'];
      else
         solver = ['QMRIDR(',int2str(s),')'];
      end
   end
   figure(fig);
   plot(log10(tot_res/norm(b)),c(kl));
   methods = strvcat(methods,solver);
   legend(methods);

end
