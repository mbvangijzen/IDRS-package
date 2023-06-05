%
% IDRS parameterised testproblems
%
% The software is distributed without any warranty.
%
% Martin van Gijzen
% Copyright (c) 2023
%

clear all;
close all;
clc;

disp('FDM discretisation of a 3D convection-diffusion-reaction problem on a unit cube');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++');

% Define system

% Defaults:
m = 60;
eps = 1e-2;
beta(1) = 0/sqrt(5);
beta(2) = 1/sqrt(5);
beta(3) = 2/sqrt(5);
r = 5;

Title = 'Parameters';
prompt = {'Points per direction:','Diffusion:','Convection x-direction:','Convection y-direction:','Convection z-direction:', 'Reaction:'};
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
b = A*sol;

disp(' ');
disp('The parameters of the problem are :');
disp(['Number of equations = ', num2str(n),';']);
disp(['Diffusion parameter = ', num2str(eps),';']);
disp(['Convection parameters = (',num2str(beta(1)),',',num2str(beta(2)),',',num2str(beta(3)),');']);
disp(['Reaction parameter = ',num2str(r),' (Note: positive reaction parameter gives negative shift to matrix);']);
disp(' ');

% Defaults for the iterative solvers:

tol = 1e-6;
maxit = 400;
M = [];

choice = 1;
scrsz = get(0,'ScreenSize');
fig = figure('Position',[scrsz(1) + scrsz(3)/2 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);
hold on;
xlabel('Number of MATVECS')
ylabel('|r|/|b|')
title(['Convergence convection-diffusion-reaction problem']);
grid on;

t = cputime;
disp('GMRES iteration...');
[x, flag, relres, iter, resvec] = gmres(A, b, [], tol, 400 );
time = cputime - t;
resvec = log10(resvec/resvec(1));
it = [0:1:length(resvec)-1];
plot(it,resvec,'k-+');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))])
disp(['Iterations: ',num2str(iter(2))]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

t = cputime;
disp('Bi-CGSTAB iteration...');
method = 2;
s = 1;
[x, flag, relres, iter, resvec] = idrs( A, b, M, s, tol, maxit, method );
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = idrs( A, b, M, s, tol, maxit );
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = idrs( A, b, M, s, tol, maxit );
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = idrs( A, b, M, s, tol, maxit  );
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = idrs( A, b, M, s, tol, maxit );
drawnow;
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = qmridr( A, b, M, s, tol, maxit );
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = qmridr( A, b, M, s, tol, maxit );
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = qmridr( A, b, M, s, tol, maxit  );
time = cputime - t;
resvec = log10(resvec/resvec(1));
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
[x, flag, relres, iter, resvec] = qmridr( A, b, M, s, tol, maxit );
drawnow;
time = cputime - t;
resvec = log10(resvec/resvec(1));
it = [0:1:length(resvec)-1];
plot(it,resvec,'g-s');
drawnow;
disp(['Final accuracy: ', num2str(norm(b-A*x)/norm(b))]);
disp(['Iterations: ',num2str(iter)]);
disp(['CPU time: ',num2str(time),'s.']);
disp(' ');

legend('GMRES', 'Bi-CGSTAB', 'IDR(1)', 'IDR(2)', 'IDR(4)', 'IDR(8)', 'QMRIDR(1)', 'QMRIDR(2)', 'QMRIDR(4)', 'QMRIDR(8)'  );
hold off;

