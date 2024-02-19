function [x, flag, relres, iter, resvec] = ...
   qmridr(A,b,s,M1,tol,maxit, in_s,in_tol,in_it,x0,omega );
%
%QMRIDR: Quasi-MinimaL Residual IDR(s) method
%   X = QMRIDR(A,B) solves the systems of linear equations A*X=B
%
%   QMRIDR is a robust implementation that imposes orthogonality conditions
%   to calculate the iteration vectors. QMRIDR is a flexible method, and
%   can be combined with IDRS as inner iterative method. In the first S iterations
%   QMRIDR is mathematically equavilent with (flexible) GMRES.
%
%   The N-by-N coefficient matrix A must be square and the right-hand
%   side column vector B must have length N.
%
%   X = QMRIDR(A,B,S) specifies the dimension of the 'shadow space'. If S = [], then
%      QMRIDR uses the default S = 4. Normally, a higher S gives faster convergence,
%      but also makes the method more expensive.
%
%   X = QMRIDR(A,B,S,M1) uses a (left) preconditioner M1. M1 must be N-by-N and invertable.
%
%   X = QMRIDR(A,B,S,M1,TOL) specifies the tolerance of the method.  If TOL is []
%      then QMRIDR uses the default, 1e-8. Note that a cheap upperbound on the residual norm
%      is used as termination criterion.
%
%   X = QMRIDR(A,B,S,M1,TOL,MAXIT) specifies the maximum number of iterations.  If
%      MAXIT is [] then QMRIDR uses the default, min(2*N,1000).
%
%   X = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S) specifies the parameter S 
%      for the inner IDRS iterations, default is S=4
%
%   X = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL) specifies the tollerance
%      for the inner IDRS iterations. Default is IN_TOL = 1e-1
%
%   X = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL,IN_IT) specifies the maximum
%      of inner IDRS iterations. Default is IN_IT = 0 (no inner iterations)
%
%   X = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL,IN_IT,X0) specifies a user defined
%      initial guess X0. Default: X0 = 0.
%
%   X = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL,IN_IT,X0,OMEGA) specifies user values for OMEGA.
%      These are used cyclicly. Default: OMEGA is empty, the omega's are computed inside
%      the algorithm using a maintaining the convergence strategy.
%
%   [X,FLAG] = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL,IN_IT,X0,OMEGA)
%      also returns an information flag:
%       FLAG = 0: required tolerance satisfied
%       FLAG = 1: no convergence to the required tolerance within maximum
%                 number of iterations
%       FLAG = 2: check RELRES, possible stagnation above required
%                 tolerance level
%       FLAG = 3: one of the iteration parameters became zero, causing break down
%
%   [X,FLAG,RELRES] = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL,IN_IT,X0,OMEGA) also
%      returns an upper bound on the relative residual norm.
%
%   [X,FLAG,RELRES,ITER] =  QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL,IN_IT,X0,OMEGA) 
%      also returns the number of iterations.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = QMRIDR(A,B,S,M1,TOL,MAXIT,IN_S,IN_TOL,IN_IT,X0,OMEGA) 
%      also returns a vector with an upper bound on the residual norm 
%      at each matrix-vector multiplication.
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

if ( nargout == 0 )
   help qmridr;
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
if nargin < 3 || isempty(s)
   s = 4;
end
if ( s > n )
   s = n;
end

% Check if preconditioner has appropriate size
if nargin < 4 || isempty(M1)
   M1 = speye(n,n);
else
   [m1,n1] = size(M1);
   if (m1 ~= n1)
      error('Preconditioner must be square.');
   end
   if (m1 ~= n)
      es = sprintf(['Preconditioner must be matrix of %d times %d.'],m,n);
      error(es);
   end
end

if nargin < 5 || isempty(tol)
   tol = 1e-6;
end

if nargin < 6 || isempty(maxit)
   maxit = min(2*n,1000);
end

if nargin < 7 || isempty(in_s)
   in_s = 4;
end

if nargin < 8 || isempty(in_tol)
   in_tol = 1e-1;
end

if nargin < 9 || isempty(in_it)
   in_it  = 0;
end

method = 1; % qmridr

if nargin < 10 || isempty(x0)
   x = zeros(n,1);
else
   if ~isequal(size(x0),[n,1])
      es = sprintf(['Initial guess must be a column vector of' ...
            ' length %d to match the problem size.'],n);
      error(es);
   end
   x = x0;
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

resvec = zeros(maxit+1,1);

% Compute initial residual:
normb = norm(b);
tolb = tol * normb;           % Relative tolerance
g = b - A*x;
normg = norm(g);
normr = normg;
resvec(1)=normr;

iter = 0;
flag = 0;
if (normg <= tolb)           % Initial guess is a good enough solution
   relres = normg/normb;
   resvec = resvec(1);
   return
end

% Define P and kappa (depending on the method)
if ( method == 1 )
   kappa = 0.7;
elseif ( method == 2 )
   kappa = 0;
   s = 1;
end
P = make_P( g, s, method );

M = zeros(s,s);
G = zeros(n,s);
g = g/normg;
W = zeros(n,s+1);
w = zeros(n,1);
om = 1;

% Last two nonzero coefficients of projected rhs:
phi_n = 0;
phi_n1 = normr;

iter = 0;
cs = zeros(s+2,1);
sn = zeros(s+2,1);

jj = 0;
flag = -1;
while ( flag < 0 )
   
   for k = 1:s+1
%
% Update counters:
      iter = iter + 1;
%
% First phase: make new vector in G-space:
      c = zeros(s+2,1); c(s+1) = 1;
      m = P_dot( P, g, s ); 
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
% Preconditioning:
      if ( in_it == 0 )
         vtilde = M1\v;
      else
         vtilde = idrs( A, v, in_s, M1, in_tol, in_it );
      end
%
% Compute new vector in space G_j
      g = A*vtilde;
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
            flag = 3
            relres = normr/normb;
            resvec = resvec(1:iter);
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
      r = zeros(s+3,1);
      r(2:s+3) = h;
%
      low = max(1,s+3-iter);
      for l = low:s+1,
%
% Apply Givens rotation.
         t = r(l);
         r(l)   =       cs(l) *t + sn(l)*r(l+1);
         r(l+1) = -conj(sn(l))*t + cs(l)*r(l+1);
      end 
%
% Form i-th rotation matrix.
      [cs(s+2) sn(s+2) r(s+2)] = rotg( r(s+2), r(s+3) );
%
% Update projected right-hand side
      phi_n    =       cs(s+2) *phi_n1;
      phi_n1   = -conj(sn(s+2))*phi_n1;
      cs(1:s+1) = cs(2:s+2); sn(1:s+1) = sn(2:s+2);
%
% Update vector:
      if ( abs(r(s+2)) < eps )
         resvec = resvec(1:iter,:);
         relres = normr/normb;
         flag = 3;
         return
      end
      w = (vtilde - W*r(1:s+1))/r(s+2);
      W(:,1:s) = W(:,2:s+1); W(:,s+1) = w;
%
% Compute solution:
      x = x + phi_n*w;
      normr = abs(phi_n1)*sqrt(jj+1);
      resvec(iter+1,:) = normr;
      if ( normr < tolb )
         flag = 0;
         break;
      elseif ( iter == maxit )
         flag = 1;
         break;
      end
   end
end

resvec = resvec(1:iter+1,:);
relres = norm(b - A*x )/normb;
if ( ( flag == 0 ) & ( normr >= tolb ) ) flag = 2; end;

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

