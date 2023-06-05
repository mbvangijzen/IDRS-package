function [x,flag,relres,iter,resvec,H]=idrs(A,b,M1,s,tol,maxit,method,x0,U0,omega)
%
%IDRS Induced Dimension Reduction method
%   X = IDRS(A,B) solves the system of linear equations A*X=B for X.  
%
%   IDRS is an efficient implementation that imposes bi-orthogonality conditions
%   to calculate the iteration vectors. It allows recyling of a previously computed
%   subspace. It is also possible to obtain spectral information about the 
%   preconditioned system matrix.
%
%   The N-by-N coefficient matrix A must be square and the right-hand
%   side column vector B must have length N.
%
%   X = IDRS(A,B,M1) uses a (left) preconditioner M1. M1 must be N-by-N and invertable.
%
%   X = IDRS(A,B,M1,S) specifies the dimension of the 'shadow space'. If S = [], then
%      IDRS uses the default S = 4. Normally, a higher S gives faster convergence,
%      but also makes the method more expensive.
%
%   X = IDRS(A,B,M1,S,TOL) specifies the tolerance of the method.  If TOL is []
%      then IDRS uses the default, 1e-8. 
%
%   X = IDRS(A,B,M1,S,TOL,MAXIT) specifies the maximum number of iterations.  If
%      MAXIT is [] then IDRS uses the default, min(2*N,1000).
%
%   X = IDRS(A,B,SHIFT,S,TOL,MAXIT,METHOD) specifies a specific method
%      if METHOD = 1 the matrix P ("shadow space") consists S of orthogonalised
%         random vectors and the "maintaining the convergence" choice for omega is used
%      if METHOD = 2 the parameters are chosen such that the method is equivalent to BICGSTAB:
%         s = 1, P equals the right-hand-side vector, and the omega's are calculated by
%         minimizing the residual norm
%      The default is method=1
%
%   X = IDRS(A,B,M1,S,TOL,MAXIT,METHOD,X0) specifies a user defined
%      initial guess X0. Default: X0 = 0.
%
%   X = IDRS(A,B,M1,S,TOL,MAXIT,METHOD,X0,U0) specifies a user defined 
%      initial search space U0. This provides a simple way to recycle information.
%      the matrix U0 should be NxS. Default: U0 = [].
%
%   X = IDRS(A,B,M1,S,TOL,MAXIT,IN_S,X0,U0,OMEGA) specifies user values for OMEGA.
%      These are used cyclicly. Default: OMEGA = [], the omega's are computed inside
%      the algorithm using a maintaining the convergence strategy.
%
%   [X,FLAG] = IDRS(A,B,M1,S,TOL,MAXIT,IN_S,X0,U0,OMEGA)
%      also returns an information flag:
%       FLAG = 0: required tolerance satisfied
%       FLAG = 1: no convergence to the required tolerance within maximum
%                 number of iterations
%       FLAG = 2: check RELRES, possible stagnation above required
%                 tolerance level
%       FLAG = 3: one of the iteration parameters became zero, causing break down
%
%   [X,FLAG,RELRES] = IDRS(A,B,M1,S,TOL,MAXIT,METHOD,X0,U0,OMEGA)
%      also returns the relative residual norm.
%
%   [X,FLAG,RELRES,ITER] = IDRS(A,B,M1,S,TOL,MAXIT,METHOD,X0,U0,OMEGA)
%      also returns the number of iterations.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = IDRS(A,B,M1,S,TOL,MAXIT,METHOD,X0,U0,OMEGA)
%      also returns a vector with the residual norm at each matrix-vector multiplication.
%
%   [X,FLAG,RELRES,ITER,RESVEC,H] = IDRS(A,B,M1,S,TOL,MAXIT,METHOD,X0,U0,OMEGA)
%      also returns a Hessenberg matrix. The eigenvalues of this matrix are Ritzvalues,
%      and converge to the eigenvalues of the preconditioned system matrix.
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

if ( nargout == 0 )
   help idrs;
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

% Check if preconditioner has appropriate size
if nargin < 3 || isempty(M1)
   M1 = speye(n,n);
else
   [m1,n1] = size(A);
   if (m1 ~= n1)
      error('Preconditioner must be square.');
   end
   if (m1 ~= n)
      es = sprintf(['Preconditioner must be matrix of %d times %d.'],m,n);
      error(es);
   end
end

% Assign default values to unspecified parameters
if nargin < 4 || isempty(s)
   s = 4;
end
if ( s > n )
   s = n;
end
if nargin < 5 || isempty(tol)
   tol = 1e-6;
end
if nargin < 6 || isempty(maxit)
   maxit = min(2*n,1000);
end

if nargin < 7 || isempty(method)
   method = 1;
end
   
if nargin < 8 || isempty(x0)
   x = zeros(n,1);
else
   if ~isequal(size(x0),[n,1])
      es = sprintf(['Initial guess must be a column vector of' ...
            ' length %d to match the problem size.'],n);
      error(es);
   end
   x = x0;
end

if nargin < 9 || isempty(U0)
   inispace = 0;
   U = zeros(n,s);
else
   if ~isequal(size(U0),[n,s])
      es = sprintf(['Initial search space must be a matrix of' ...
            ' %d times %d to match the problem size.'],n,s);
      error(es);
   end
   U = U0;
   inispace = 1;
end

user_omega = 1;
if ( nargin < 10 || isempty(omega)  ) 
   user_omega = 0;
   omega = [];
end

if nargin > 10
   es = sprintf(['Too many input parameters']);
   error(es);
end

out_H = ( nargout == 6 );
resvec = zeros(maxit+1,1);

% Compute initial residual:
normb = norm(b);
tolb = tol * normb;           % Relative tolerance
r = b - A*x;
normr = norm(r);
resvec(1)=normr;

iter = 0;                 
flag = 0;
if (normr <= tolb)           % Initial guess is a good enough solution
   relres = normr/normb;
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
P = make_P( r, s, method );

G = zeros(n,s); M = zeros(s,s); 
om = 1;
iter = 0;
jj = 0;
ii = 0;

alpha = zeros(s,1);
beta  = zeros(s,1);
gamma = zeros(s,1);
f     = zeros(s,1);

if (out_H) H = zeros(maxit+1,maxit); end;

% Main iteration loop, build G-spaces:
while ( normr > tolb && iter < maxit )  

% New righ-hand size for small system:
   f = P_dot(P,r,s);
   for k = 1:s 
      ii = ii + 1;

% Compute new U(:,k) and G(:,k), G(:,k) is in space G_j
      v = r;
      if ( jj > 0 ) 

% Solve small system and make v orthogonal to P:
         gamma = M(k:s,k:s)\f(k:s); 
         v = v - G(:,k:s)*gamma;
         U(:,k) = U(:,k:s)*gamma + om*(M1\v);
         if ( out_H ) H(ii-s:ii-k,ii)   =  -gamma./beta(k:s); end;
%
      elseif ~( inispace )
         U(:,k) = M1\v;
      end
      G(:,k) = A*U(:,k);
%
% Bi-Orthogonalise the new basis vectors:
      mu = P_dot(P,G(:,k),s);
      if ( k > 1 )
         alpha = M(1:k-1,1:k-1)\mu(1:k-1);
         G(:,k) = G(:,k) - G(:,1:k-1)*alpha;
         U(:,k) = U(:,k) - U(:,1:k-1)*alpha;
         mu(k:s) = mu(k:s) - M(k:s,1:k-1)*alpha;
      end
      M(k:s,k) = mu(k:s);
%
      if ( k > 1 && out_H ) H(ii-k+1:ii-1,ii) = alpha(1:k-1)./beta(1:k-1); end;

% Break down?
      if ( M(k,k) == 0 )
         flag = 3;
         relres = normr/normb;
         resvec = resvec(1:iter);
         return;
      end
%
% Make r orthogonal to p_i, i = 1..k 
      beta(k) = f(k)/M(k,k);
      r = r - beta(k)*G(:,k);
      x = x + beta(k)*U(:,k);
%
% New f = P'*r (first k  components are zero)
      if ( k < s ) 
         f(k+1:s)   = f(k+1:s) - beta(k)*M(k+1:s,k);
      end
%
      if ( out_H )
         H(ii,ii) = 1/beta(k);
         l = max(1,ii-s);
         H(l+1:ii+1,ii) = H(l+1:ii+1,ii) - H(l:ii,ii);
         H(l:ii+1,ii)   = H(l:ii+1,ii)/om;
      end
%
      iter = iter + 1;
      normr = norm(r);
      resvec(iter+1) = normr;
%      
      if ( normr < tolb | iter == maxit ) 
         break
      end
   end 
%
   if ( normr < tolb | iter == maxit )
      break
   end
%
% Now we have sufficient vectors in G_j to compute residual in G_j+1
% Note: r is already perpendicular to P so v = r
%
   jj = jj+1;

% Preconditioning:
   v = M1\r;
   t = A*v;

% Computation of a new omega
   if ( user_omega )
      ind = rem(jj,length(omega));
      if ( ind == 0 ) ind = length(omega); end;
      om = omega(ind);
   else
      om = comp_om( t, r, kappa );
   end

   if ( om == 0 )
      flag = 3;
      relres = normr/normb;
      resvec = resvec(1:iter);
      return;
   end
%
   r = r - om*t;
   x = x + om*v;
   normr = norm(r);
%
   iter = iter + 1;
   resvec(iter+1) = normr;
   
end; %while

resvec = resvec(1:iter+1);
if ( out_H ) H = H(1:ii+1,1:ii); end;
relres = norm(b - A*x)/normb;
if ( relres < tol ) 
   flag = 0;
elseif ( iter == maxit )
   flag = 1;
else
   flag = 2;
end

return

