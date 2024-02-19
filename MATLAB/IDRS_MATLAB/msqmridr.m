function [x, flag, relres, iter, resvec]=msqmridr(A,b,sigma,s,M1,tol,maxit,in_s,in_tol,in_it,omega);
%
%MSQMRIDR Multi Shift Induced Dimension Reduction method
%   X = MSQMRIDR(A,B,SHIFT) simultaneously solves the sequence of shifted systems 
%             (A-SHIFT(i) I) X_i = B or (A-SHIFT(i) M)X_i = B for X_i.
%
%   MSQMRIDR is a robust implementation that imposes orthogonality conditions
%   to calculate the iteration vectors. MSQMRIDR is a flexible method, and
%   can be combined with MSIDRS as inner iterative method. 
%
%   The N-by-N coefficient matrix A must be square and the right-hand
%   side column vector B must have length N.
%
%   X = MSQMRIDR(A,B,SHIFT) specifies the shifts. If SHIFT = [], no shift is applied,
%      which had the same  effect as SHIFT = 0.
%
%   X = MSQMRIDR(A,B,SHIFT,S) specifies the dimension of the 'shadow space'. If S = [], then
%      MSQMRIDR uses the default S = 4. Normally, a higher S gives faster convergence,
%      but also makes the method more expensive.
%
%   X = MSQMRIDR(A,B,SHIFT,S,M) specifies the mass matrix.  If M is []
%      then MSQMRIDR uses M = I.
%
%   X = MSQMRIDR(A,B,SHIFT,S,M,TOL) specifies the tolerance of the method.  If TOL is []
%      then MSQMRIDR uses the default, 1e-6. Note that a cheap upperbound on the residual norms
%      is used as termination criterion.
%
%   X = MSQMRIDR(A,B,SHIFT,S,M,TOL,MAXIT) specifies the maximum number of iterations.  If
%      MAXIT is [] then MSQMRIDR uses the default, min(2*N,1000).
%
%   X = MSQMRIDR(A,B,SHIFT,S,M,TOL,MAXIT,IN_S) specifies the parameter S
%      for the inner MSIDRS iterations, default is S=4
%
%   X = MSQMRIDR(A,B,SHIFT,S,M,TOL,MAXIT,IN_S,IN_TOL) specifies the tollerance
%      for the inner MSIDRS iterations. Default is IN_TOL = 1e-1
%
%   X = MSQMRIDR(A,B,SHIFT,S,TOL,M,MAXIT,IN_S,IN_TOL,IN_IT) specifies the maximum
%      of inner MSIDRS iterations. Default is IN_IT = 0 (no inner iterations)
%
%   X = MSQMRIDR(A,B,SHIFT,S,TOL,MAXIT,IN_S,IN_TOL,IN_IT,OMEGA) specifies user values for OMEGA.
%      These are used cyclicly. Default: OMEGA is empty, the omega's are computed inside
%      the algorithm using a maintaining the convergence strategy.
%
%   [X,FLAG] = MSQMRIDR(A,B,SHIFT,S,TOL,MAXIT,IN_S,IN_TOL,IN_IT,OMEGA)
%      also returns an information flag:
%       FLAG = 0: required tolerance satisfied
%       FLAG = 1: no convergence to the required tolerance within maximum
%                 number of iterations
%       FLAG = 2: check RELRES, possible stagnation above required
%                 tolerance level
%       FLAG = 3: one of the iteration parameters became zero, causing break down
%
%   [X,FLAG,RELRES] = MSQMRIDR(A,B,SHIFT,S,TOL,MAXIT,IN_S,IN_TOL,IN_IT,OMEGA) also
%      returns an upper bound on the relative residual norms.
%
%   [X,FLAG,RELRES,ITER] = MSQMRIDR(A,B,SHIFT,S,TOL,MAXIT,IN_S,IN_TOL,IN_IT,OMEGA)
%      also returns the number of iterations.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = MSQMRIDR(A,B,SHIFT,S,TOL,MAXIT,IN_S,IN_TOL,IN_IT,OMEGA)
%      also returns a vector with an upper bound on the residual norm
%      at each matrix-vector multiplication.
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%


if ( nargout == 0 )
   help msqmridr;
   return
end

% Check for an acceptable number of input arguments
if nargin < 2
   error('Not enough input arguments.');
end

% Check matrix and right hand side vector inputs have appropriate sizes
[m,n] = size(A);
if (m ~= n)
   error('Matrix must be square.');
end
if ~isequal(size(b),[m,1])
   es = sprintf(['Right hand side must be a column vector of' ...
         ' length %d to match the coefficient matrix.'],m);
   error(es);
end

% Assign default values to unspecified parameters
if nargin < 3 || isempty(sigma)
   sigma = 0;
   n_sigma = 1;
else
   n_sigma = length(sigma);
end

if nargin < 4 || isempty(s)
   s = 4;
end
if ( s > n )
   s = n;
end

if nargin < 5 || isempty(M1)
   mass = 0;
else
   mass = 1;
end

if nargin < 6 || isempty(tol)
   tol = 1e-6;
end

if nargin < 7 || isempty(maxit)
   maxit = min(2*n,1000);
end

method = 1; % qmridr

if nargin < 8 || isempty(in_s)
   in_s = 4;
end
if nargin < 9 || isempty(in_tol)
   in_tol = 1e-1;
end
if nargin < 10 || isempty(in_it)
   in_it  = 0;
end

user_omega = 1;
if ( nargin < 11 || isempty(omega)  )
   user_omega = 0;
   omega = [];
end

if nargin > 11
   es = sprintf(['Too many input parameters']);
   error(es);
end

% Compute initial residual:
normb = norm(b);
tolb = tol * normb;           % Relative tolerance

x = zeros(n,n_sigma);
% Check for zero rhs:
iter = 0;
flag = 0;
if (normb == 0)              % Solution is nulvector
   resvec = zeros(1,n_sigma);
   relres = zeros(1,n_sigma);
   return
end

% Define P and kappa (depending on the method)
if ( method == 1 )
   kappa = 0.7;
elseif ( method == 2 )
   kappa = 0;
   s = 1;
end
P = make_P( b, s, method );

% Initialize output paramater relres
relres = ones(1,n_sigma);
normr = normb*ones(1,n_sigma);
resvec = zeros(maxit,n_sigma);
resvec(1,:) = normr;

g = b/normb;
G = zeros(n,s); M = zeros(s,s);

iter = 0;

W  = zeros(n,s+1,n_sigma);
w  = zeros(n,1);
om = sqrt(2);

% Last two nonzero coefficients of projected rhs:
phi_n =  zeros(n_sigma);
phi_n1 = zeros(n_sigma);
phi_n1 = normr;

cs = zeros(s+2,n_sigma);
sn = zeros(s+2,n_sigma);

converged = zeros(n_sigma,1);

% Main iteration loop, build G-spaces:
jj = 0;

while ( max(normr) > tolb && iter < maxit )

   for k = 1:s+1
%
% Update counters:
      iter = iter + 1;
%
% First phase: make new vector in G-space:
      c = zeros(s+2,1); c(s+1) = 1;
      m = P_dot(P,g,s); 
      if iter > s
%
% Construct v orthogonal to P
         gamma = M\m;
         v = g - G*gamma;
	 c(1:s) = -gamma;
      else  
%
% First s steps: Arnoldi
         v = g;
      end
      M(:,1:s-1) = M(:,2:s); M(:,s) = m;
      G(:,1:s-1) = G(:,2:s); G(:,s) = g; 
%
% Inner iterations?
      if ( in_it == 0 )
         vtilde = v;
      else
         sigma_0(1:n_sigma) = sigma;
         sigma_0(n_sigma+1) = 0.;
         if ( mass ) 
            [V_tilde, in_flag, in_relres, in_iter, in_resvec, cf]  = msidrs( A, v, sigma_0, in_s, M1, in_tol, in_it );
         else
            [V_tilde, in_flag, in_relres, in_iter, in_resvec, cf]  = msidrs( A, v, sigma_0, in_s, [],  in_tol, in_it );
         end 
         eta = cf(1:n_sigma)/cf(n_sigma+1);
         vtilde = V_tilde(:,n_sigma+1);
      end
%
% Compute new vector in space G_j
      if ( mass ) 
         g = A*(M1\vtilde);
      else
         g = A*vtilde;
      end
%
% New G-space? 
      if ( k == s+1 )
% Computation of a new omega
         if ( user_omega )
            ind = rem(jj,length(omega));
            if ( ind == 0 ) ind = length(omega); end;
            om = omega(ind);
         else
            om =  comp_om( g, v, kappa );
         end
         if ( abs(om) < eps )
            info = 3
            return
         end
         jj = jj+1;
      end
      g = g - v/om;
%
% Orthogonalisation (2 times classical Gram-Schmidt)
      h = c/om;
      if ( k < s+1 )
         beta  = G(:,s-k+1:s)'*g; g = g - G(:,s-k+1:s)* beta;
         cbeta = G(:,s-k+1:s)'*g; g = g - G(:,s-k+1:s)*cbeta;
         beta = beta + cbeta;
         h(s+1-k+1:s+1) = h(s+1-k+1:s+1)+beta;
      end
%
% Normalise
      normg = norm(g);
      g = g/normg;
      h(s+2) = normg;
%
% Store the information
      for k_sigma = 1:n_sigma
         if ( ~converged(k_sigma) )
            r_sigma = zeros(s+3,1);
            if ( in_it == 0 )
               r_sigma(2:s+3) = h - sigma(k_sigma)*c;
            else
               r_sigma(2:s+3) = eta(k_sigma)*h - (eta(k_sigma)-1.)*c;
            end
%
            low = max(1,s+3-iter);
            for l = low:s+1
%
% Apply Givens rotation.
               t = r_sigma(l); 
               r_sigma(l)   =       cs(l,k_sigma) *t + sn(l,k_sigma)*r_sigma(l+1);
               r_sigma(l+1) = -conj(sn(l,k_sigma))*t + cs(l,k_sigma)*r_sigma(l+1);
            end
%
% Form i-th rotation matrix.
            [cs(s+2,k_sigma) sn(s+2,k_sigma) r_sigma(s+2)] = rotg( r_sigma(s+2), r_sigma(s+3) );
%
% Update projected right-hand side
            phi_n(k_sigma)    =       cs(s+2,k_sigma) *phi_n1(k_sigma);
            phi_n1(k_sigma)   = -conj(sn(s+2,k_sigma))*phi_n1(k_sigma);
            cs(1:s+1,k_sigma) = cs(2:s+2,k_sigma); sn(1:s+1,k_sigma) = sn(2:s+2,k_sigma);
%
% Update vector:
            if ( abs(r_sigma(s+2)) < eps )
               resvec = resvec(1:iter,:);
               flag = 3;
               return
            end
            if ( in_it > 0 ) vtilde = V_tilde(:,k_sigma); end;
            w = (vtilde - W(:,:,k_sigma)*r_sigma(1:s+1))/r_sigma(s+2);
            W(:,1:s,k_sigma) = W(:,2:s+1,k_sigma); W(:,s+1,k_sigma) = w;
%
% Compute solution:
            x(:,k_sigma) = x(:,k_sigma) + phi_n(k_sigma)*w;
            normr(k_sigma) = abs(phi_n1(k_sigma))*sqrt(jj+1);
	    converged(k_sigma) = normr(k_sigma) < tolb;
            resvec(iter+1,k_sigma) = normr(k_sigma);
         end 
      end
      if ( max(normr) < tolb | iter == maxit )
         break
      end
   end
end

resvec = resvec(1:iter+1,:);
if ( mass )
   for k_sigma = 1:n_sigma
      relres(k_sigma) = norm(b - A*(M1\x(:,k_sigma)) + sigma(k_sigma)*x(:,k_sigma))/normb;
      x(:,k_sigma) = M1\x(:,k_sigma);
   end
else
   for k_sigma = 1:n_sigma
      relres(k_sigma) = norm(b - A*x(:,k_sigma) + sigma(k_sigma)*x(:,k_sigma))/normb;
   end
end 

if ( max(relres) < tol )
   flag = 0;
elseif ( iter == maxit )
   flag = 1;
else
   flag = 2;
end

return

function [c,s,r] = rotg(a,b);
% 
% ROTG: construct and apply Givens rotation
% ROTG returns c,s,r such that
%        | c       s || x |   | r |
%        |           ||   | = |   |
%        |-conj(s) c || y | = | 0 |
%
% ROTG is a translation of the BLAS routine CROTG
%
if abs(a) < eps
   c = 0;
   s = 1;
   r = b;
else
   scale = abs(a) + abs(b);
   rho = scale*sqrt(abs(a/scale)^2 +  abs(b/scale)^2);
   alpha = a/abs(a);
   c = abs(a)/rho;
   s = alpha*conj(b)/rho;
   r = alpha*rho;
end
%
return
