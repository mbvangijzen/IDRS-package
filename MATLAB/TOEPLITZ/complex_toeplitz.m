%
% 1D convection diffusion test problem
%
% The software is distributed without any warranty.
%
% Martin van Gijzen
% Copyright (c) February 2014
%

clear all;
close all;
clc;

disp('Complex Toeplitz  problem');
disp('+++++++++++++++++++++++++');

% Define system

% Defaults:
N = 200;
e = ones(N,1);
im = sqrt(-1);
c(1) = 3.6*im;
c(2) = 4.;
c(3) = 1.;
c(4) = 0.7;

A = spdiags([c(1)*e c(2)*e c(3)*e c(4)*e], [-1,0,2,3], N, N );

b = ones(N,1)*im;

disp('Accuracy test problem, see TOMS paper (2011), example 2.');
disp('The termination criterion is ||r||/||b|| < 1e-12.');
disp(' ');

% Parameters for the iterative solvers:

tol = 1e-12;
maxit = 1000;

choice = 1;
scrsz = get(0,'ScreenSize');
fig = figure('Position',[scrsz(1) + scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
hold on;
xlabel('Number of MATVECS')
ylabel('|r|/|b|')
title(['IDR accuracy test']);
grid on;

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

s = 16;
t = cputime;
disp('IDR(16) iteration...');
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

s = 32;
t = cputime;
disp('IDR(32) iteration...');
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

s = 64;
t = cputime;
disp('IDR(64) iteration...');
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

s = 1;
t = cputime;
disp('QMRIDR(1) iteration...');
% Compute solution
[x, flag, relres, iter, resvec] = qmridr( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-+');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 2;
t = cputime;
disp('QMRIDR(2) iteration...');
% Initial iterations to compute ritz values:
[x, flag, relres, iter, resvec] = qmridr( A, b, s, [], tol, maxit );
time = cputime - t;
figure(fig);
resvec = log10(resvec/norm(b));
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-x');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 4;
t = cputime;
disp('QMRIDR(4) iteration...');
% Compute solution
[x, flag, relres, iter, resvec] = qmridr( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-*');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 8;
t = cputime;
disp('QMRIDR(8) iteration...');
% Compute solution
randn('seed',0);
%x0 = randn(n,1);
[x, flag, relres, iter, resvec] = qmridr( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-s');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 16;
t = cputime;
disp('QMRIDR(16) iteration...');
% Compute solution
randn('seed',0);
%x0 = randn(n,1);
[x, flag, relres, iter, resvec] = qmridr( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-s');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 32;
t = cputime;
disp('QMRIDR(32) iteration...');
% Compute solution
randn('seed',0);
%x0 = randn(n,1);
[x, flag, relres, iter, resvec] = qmridr( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-s');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

s = 64;
t = cputime;
disp('QMRIDR(64) iteration...');
% Compute solution
randn('seed',0);
%x0 = randn(n,1);
[x, flag, relres, iter, resvec] = qmridr( A, b, s, [], tol, maxit );
time = cputime - t;
resvec = log10(resvec/norm(b));
figure(fig);
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-s');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

legend('BiCGSTAB', 'IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)', 'IDR(16)', 'IDR(32)', 'IDR(64)', ...
'QMRIDR(1)', 'QMRIDR(2)', 'QMRIDR(4)', 'QMRIDR(8)', 'QMRIDR(16)', 'QMRIDR(32)', 'QMRIDR(64)');
hold off;

