function [x, flag, relres, iter, resvec,colfac]=msidrs(A,b,sigma,s,M1,tol,maxit,method,omega );
%
%MSIDRS: Multi-shift IDR(s) method
%   X = MSIDRS(A,B) simultaneously solves the systems of linear equations 
%                  (A-SHIFT(i)I)*X=B, i= 1, 2, ...,  for X.
%
%   MSIDRS is an efficient implementation that imposes bi-orthogonality conditions
%   to calculate the iteration vectors. The residuals of the shifted problems
%   are collinear. This property can be exploited by MSQMRIDR by using MSIDRS as inner
%   iterative method.
%
%   The N-by-N coefficient matrix A must be square and the right-hand
%   side column vector B must have length N. 
%
%   X = MSIDRS(A,B,SHIFT) specifies the shifts. If SHIFT = [], no shift is applied,
%   which had the same  effect as SHIFT = 0;   
%
%   X = MSIDRS(A,B,SHIFT,S) specifies the dimension of the 'shadow space'. If S = [], then
%   MSIDRS uses the default S = 4. Normally, a higher S gives faster convergence,
%   but also makes the method more expensive.
%
%   X = MSIDRS(A,B,SHIFT,S,M) specifies the mass matrix.  If M is []
%   then MSIDRS uses M = I.
%
%   X = MSIDRS(A,B,SHIFT,S,M,TOL) specifies the tolerance of the method.  If TOL is []
%   then MSIDRS uses the default, 1e-6.
%
%   X = MSIDRS(A,B,SHIFT,S,M,TOL,MAXIT) specifies the maximum number of iterations.  If
%   MAXIT is [] then MSIDRS uses the default, min(2*N,1000).
%
%   X = MSIDRS(A,B,SHIFT,S,M,TOL,MAXIT,METHOD) specifies a specific method
%      if METHOD = 1 the matrix P ("shadow space") consists S of orthogonalised 
%         random vectors and the "maintaining the convergence" choice for omega is used
%      if METHOD = 2 the bicgstab parameters are used: 
%         s = 1 and P equals the right-hand-side vector, and the omega's are calculated by
%         minimizing the resdual norm
%      The default is method=1
%
%   X = MSIDRS(A,B,SHIFT,S,M,TOL,MAXIT,METHOD,OMEGA) defines a sequence of
%      user-defined parameters OMEGAs. The OMEGAs are used cyclicly.
%      By default the "maintaining the convergence" method is used for OMEGA.
%
%   [X,FLAG] = MSIDRS(A,B,SHIFT,S,M,TOL,MAXIT,METHOD,OMEGA)
%   also returns an information flag:
%       FLAG = 0: required tolerance satisfied
%       FLAG = 1: no convergence to the required tolerance within maximum
%                 number of iterations
%       FLAG = 2: check RELRES, possible stagnation above required
%                 tolerance level
%       FLAG = 3: one of the iteration parameters became zero, causing break down
%
%   [X,FLAG,RELRES] = MSIDRS(A,B,S,M,TOL,MAXIT,METHOD,OMEGA) also
%      returns the relative residual norms.
%
%   [X,FLAG,RELRES,ITER] = MSRIDRS(A,B,SHIFT,S,M,TOL,MAXIT,METHOD,OMEGA) also returns
%      the number of iterations.
%
%   [X,FLAG,RELRES,ITER,RESVEC] = MSIDRS(A,B,SHIFT,S,M,TOL,MAXIT,METHOD,OMEGA) also returns
%      a vector with, for every shift, the residual norm at each matrix-vector multiplication.
%
%   [X,FLAG,RELRES,ITER,RESVEC,COLFAC] = MSIDRS(A,B,SHIFT,S,M,TOL,MAXIT,METHOD,OMEGA)
%      also returns the colinearity factors (for use in MSQMRIDR)
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

if ( nargout == 0 )
   help msidrs;
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

if nargin < 8 || isempty(method)
   method = 1;
end

user_omega = 1;
if ( nargin < 9 || isempty(omega)  )
   user_omega = 0;
   omega = [];
end

if nargin > 9
   es = sprintf(['Too many input parameters']);
   error(es);
end

out_colfac = 0;
if ( nargout == 6 )
   out_colfac = 1;
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

g = b;
G = zeros(n,s); M = eye(s,s);

iter = 0;

W = zeros(n,s,n_sigma);
w  = zeros(n,1);
v0  = zeros(n,1);
om = sqrt(2);
om_sigma = om./(1 - om*sigma);

phi =  ones(1,n_sigma);

l_sigma = zeros(s+2,n_sigma);

% Main iteration loop, build G-spaces:
jj = 0;

while ( max(normr) > tolb & iter < maxit )  

%
% Compute intermediate vectors:
%
% New right-hand size for small system:
   f = P_dot(P,g,s);
   for k = 1:s 
      iter = iter + 1;
%
% Solve small system and make v orthogonal to P:
      v = g;
      if ( jj > 0 )
         gamma      = zeros(s,1);
         gamma(k:s) = M(k:s,k:s)\f(k:s); 
         v = v - G(:,k:s)*gamma(k:s);
      end
      M(k:s,k) = f(k:s);
      G(:,k) = g;
%
% Compute new g = (I - om A) v 
      if ( mass ) 
         t = A*(M1\v);
      else
         t = A*v;
      end
      g = v - om*t;
%
% Bi-Orthogonalise the new basis vectors: 
      f = P_dot(P,g,s);
      alpha      = zeros(s,1);
      alpha(1:k) = M(1:k,1:k)\f(1:k);
      g = g - G(:,1:k)*alpha(1:k);
      f(k+1:s) = f(k+1:s) - M(k+1:s,1:k)*alpha(1:k);
      normg = norm(g);
%
% Now form the next column of the Hessenberg matrix
      h = zeros(s+3,1);
      if ( jj > 0 ) h(1:s-k+1) = -gamma(k:s); end;
      h(s+2-k+1:s+2) = h(s+2-k+1:s+2)-alpha(1:k);
      h(s+2) = h(s+2)+1.;
      h(s+3) = -1;
      h = h/om;
      for k_sigma = 1:n_sigma
         u_sigma = h;
         if ( jj > 0 ) 
            u_sigma(1:s-k+1) = u_sigma(1:s-k+1) + sigma(k_sigma)*gamma(k:s);
         end;
         u_sigma(s+2) = u_sigma(s+2)-sigma(k_sigma);

         l_sigma(1:s+1,k_sigma) = l_sigma(2:s+2,k_sigma);
         for j = 2:s+2,
            u_sigma(j) = u_sigma(j) - l_sigma(j-1,k_sigma)*u_sigma(j-1);
         end
%
% Next pivot, update solution of projected system
         l_sigma(s+2,k_sigma) = -1/(om*u_sigma(s+2));

%
% New Update vector:
% Corresponds with Dimension Reduction step
         w = v - om_sigma(k_sigma)*v0*u_sigma(s+2-k);
% Corresponds with space j
         if ( k > 1 ) w = w - W(:,1:k-1,k_sigma)*u_sigma(s+2-k+1:s+1); end;
% Corresponds with space j+1
         w = w - W(:,k:s,k_sigma)*u_sigma(1:s+2-k-1); 
% Scale and store
         w = w/u_sigma(s+2);
         W(:,k,k_sigma) = w;
%
% Compute solution:
         x(:,k_sigma) = x(:,k_sigma) + phi(k_sigma)*w;
% Update rhs of projected system:
         phi(k_sigma)    =  -l_sigma(s+2,k_sigma)*phi(k_sigma);
% Compute residual norms
         normr(k_sigma) = normg*abs(phi(k_sigma));
      end 
      resvec(iter+1,:) = normr;
      if ( max(normr) < tolb | iter == maxit )
         break
      end
%
   end
   if ( max(normr) < tolb | iter == maxit )
      break
   end
%
% Dimension reduction step
   iter = iter + 1;
   jj = jj +1;
   v = g;
%
   if ( mass ) 
      t = A*(M1\v); 
   else
      t = A*v;
   end

% Computation of a new omega
   if ( user_omega )
      ind = rem(jj,length(omega));
      if ( ind == 0 ) ind = length(omega); end;
      om = omega(ind);
   else
      om = comp_om( t, g, kappa );
   end
   g = g - om*t;
   normg = norm(g);
%
%
% It is not necessary to form a new colum of the Hessenberg matrix, it is simply
%        h_sigma(s+2) = 1/om
%        h_sigma(s+3) = -1/om;
% For u_sigma we therefore get
%        u_sigma(s+2) = 1/om - sigma;
%
   for k_sigma = 1:n_sigma
% Compute pivots
      l_sigma(1:s+1,k_sigma) = l_sigma(2:s+2,k_sigma); 
      l_sigma(s+2,k_sigma)   = -1/(1-om*sigma(k_sigma));

% Compute u_sigma(s+2), these are omega's for the shifted systems:
      om_sigma(k_sigma) = om/(1-om*sigma(k_sigma)); % This is u_sigma(s+2)
%
% Update vector. Note that this one is special: 
% It is collinear for all the shifts, with collinearity factor om_sigma
      v0 = v;
%
% Compute solution:
      x(:,k_sigma) = x(:,k_sigma) + (phi(k_sigma)*om_sigma(k_sigma))*v0;
% Update rhs of projected system:
      phi(k_sigma)      = - l_sigma(s+2,k_sigma)*phi(k_sigma);
% Compute residual norms
      normr(k_sigma) = normg*abs(phi(k_sigma));
      resvec(iter+1,k_sigma) = normr(k_sigma);
   end

end; %while

resvec = resvec(1:iter+1,:);
if ( mass )
   for k_sigma = 1:n_sigma
      relres(k_sigma) = norm(b - A*(M1\x(:,k_sigma)) + sigma(k_sigma)*x(:,k_sigma))/normb;
   end
else
   for k_sigma = 1:n_sigma
      relres(k_sigma) = norm(b - A*x(:,k_sigma) + sigma(k_sigma)*x(:,k_sigma))/normb;
   end
end

% Scale back the solution (not if msidrs is used as inner iterative method):
if ( mass & out_colfac == 0 )
   for k_sigma = 1:n_sigma
      x(:,k_sigma) = M1\x(:,k_sigma);
   end
end

if ( max(relres) < tol )
   flag = 0;
elseif ( iter == maxit )
   flag = 1;
else
   flag = 2;
end
colfac = phi;

return

