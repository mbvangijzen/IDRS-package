%
% IDRS parameterised testproblems
%
% The software is distributed without any warranty.
%
% Martin van Gijzen
% Copyright (c) 2024
%

clear all;
close all;
clc;

disp('FDM discretisation of sound propagation in an acoustically dead room');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');

% Define system
im = sqrt(-1);
grid = 2;
m = 16*2^(grid-1);
L = 4;
sound_velocity = 340;
freq = 1e2*2^(grid-1);
k = 2*pi*freq/sound_velocity;

% Generate matrix
h = L/(m-1);

n = m*m*m;
S = gallery('tridiag',m,-1/h^2,2/h^2,-1/h^2);

% Correct for boundary condition:
S(1,2) = -2/h^2;
S(1,1) = S(1,1) + 2*k*im/h;
S(m-1,m) = -2/h^2;
S(m,m) = S(m,m) + 2*k*im/h;
Is = speye(m,m);
I = speye(n,n);
A = kron(kron(Is,Is),S) + kron(kron(Is,S),Is)+ kron(kron(S,Is),Is) -k^2*I;

l = m/2 +(m/2-1)*m +(m/2-1)*m^2;
b = zeros(n,1);
b(l) = 1/(h^3);

% Defaults for the iterative solvers:

tol = 1e-6;
maxit = 1000;

choice = 1;
scrsz = get(0,'ScreenSize');
fig = figure('Position',[scrsz(1) + scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
hold on;
xlabel('Number of MATVECS')
ylabel('|r|/|b|')
title('Acoustically dead room problem');
grid on;

t = cputime;
disp('GMRES iteration...');
[x, flag, relres, iter, resvec] = gmres(A, b, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'k-+');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))])
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

t = cputime;
disp('BiCGSTAB iteration...');
[x, flag, relres, iter, resvec] = bicgstab(A, b, tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'b-+');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))])
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 1;
t = cputime;
disp('IDR(1) iteration...');
% Compute solution
[x, flag, relres, iter, resvec] = idrs( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'r-+');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 2;
t = cputime;
disp('IDR(2) iteration...');
% Initial iterations to compute ritz values:
[x, flag, relres, iter, resvec] = idrs( A, b, s, [], tol, maxit );
time = cputime - t;
figure(fig);
resvec = log10(resvec/norm(b));
it = [0:1:length(resvec)-1];
plot(it,resvec,'r-x');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 4;
t = cputime;
disp('IDR(4) iteration...');
% Compute solution
[x, flag, relres, iter, resvec] = idrs( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'r-*');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 8;
t = cputime;
disp('IDR(8) iteration...');
% Compute solution
randn('seed',0);
%x0 = randn(n,1);
[x, flag, relres, iter, resvec] = idrs( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'r-s');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

legend('GMRES', 'BiCGSTAB', 'IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)' );
hold off;

x = reshape(x,m,m,m);
y = squeeze(real(x(m/2,:,:)));

figure;
surf(y);
