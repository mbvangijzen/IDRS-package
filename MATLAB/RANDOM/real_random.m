%
% This programme illustrates the use of the idrs-routines on an academic real dense system.
% The rhs-vector is such that the solution is equal to one.
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
% 


clear all;
close all;
clc;

disp('Real random test problem');
disp('++++++++++++++++++++++++');

% Define system

% Defaults:
N = 100;

rng('default');
A = rand(N,N);

b = zeros(N,1);
x = ones(N,1);
b = A*x;

disp('Real random problem, N=100 gridpoints.');
disp('The termination criterion is ||r||/||b|| < 1e-6.');
disp(' ');

% Parameters for the iterative solvers:

tol = 1e-6;
maxit = 1000;

choice = 1;
scrsz = get(0,'ScreenSize');
fig = figure('Position',[scrsz(1) + scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
hold on;
xlabel('Number of MATVECS')
ylabel('|r|/|b|')
title(['IDR real random matrix']);
grid on;

t = cputime;
disp('GMRES iteration...');
[x, flag, relres, iter, resvec] = gmres(A, b, [], tol, N );
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

legend('GMRES', 'BiCGSTAB', 'IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)', 'QMRIDR(1)', 'QMRIDR(2)', 'QMRIDR(4)', 'QMRIDR(8)');
hold off;

