%
% Multishift acousticly dead room test problem
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2024 Martin van Gijzen
%

clear all;
close all;
clc;

disp('FDM discretisation of sound propagation in an acoustically dead room');
disp('Note that this example breaks down if real random P is used as shadow vectors');
disp('Change the routine make_P in the IDRS_MATLAB directory to solve this problem');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');

% Define system
im = sqrt(-1);
grid = 2;
m = 16*2^(grid-1);
L = 4;
sound_velocity = 340;
freq = 1e2*2^(grid-1);

% Generate matrix
h = L/(m-1);

n = m*m*m;
T = gallery('tridiag',m,-1/h^2,2/h^2,-1/h^2);

% Correct for boundary condition:
T(1,2) = -2/h^2;
T(m-1,m) = -2/h^2;
Is = speye(m,m);
K = kron(kron(Is,Is),T) + kron(kron(Is,T),Is)+ kron(kron(T,Is),Is);

D = zeros(m,1);
D(1) = 2/h;
D(m) = 2/h;
D = spdiags(D,0,m,m);
C = kron(kron(Is,Is),D) + kron(kron(Is,D),Is)+ kron(kron(D,Is),Is); 

I = speye(n,n);
O = sparse(n,n);

A = [-C, -K; I, O];
l = m/2 +(m/2-1)*m +(m/2-1)*m^2;
b = zeros(2*n,1);
b(l) = -1/(h^3);

% Defaults for the iterative solvers:

tol = 1e-6;
maxit = 1000;

Title = 'Parameters';
variant = 'idrs';
s = '4';
tol = '1e-6';
maxit = '10000';
in_s = '4';
in_tol = '1e-1';
in_it = '0';

contin = 1;

freq = [100, 200];
shift = (2*pi*im/sound_velocity).*freq;

Title = 'Parameters';
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
   prompt = { 'Method: msidrs/msbicgstab/msqmridr',...
              'Parameter s:', ...
              'Tolerance:',...
              'Maximum iterations ',...
              'Inner s:',...
              'Inner tolerance:',...
              'Inner iterations:'};
   defaults = { variant, s, tol, maxit, in_s, in_tol, in_it};
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

   method = 0;
   if ( strcmpi(variant,'msidrs' ) )
      method = 1;
   elseif ( strcmpi(variant,'msbicgstab' ) )
      method = 2;
      s = '1';
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

   t = cputime;
   if ( method == 1 || method == 2 )
      [x,flag,relres,iter,resvec] = ...
          msidrs(A,b,shift,str2num(s),[],str2num(tol),str2num(maxit),method);
   else
      [x,flag,relres,iter,resvec] = ...
          msqmridr(A,b,shift,str2num(s),[],str2num(tol),str2num(maxit),...
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
      y = reshape(x(n+1:2*n,i_shift),m,m,m);
      y = squeeze(real(y(m/2,:,:)));
      figure(fig2);
      surf(y);
      drawnow;
      pause(1);
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

