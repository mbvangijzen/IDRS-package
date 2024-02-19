function [x, flag, relres, iter, resvec, U0, omega] = idrs_ritz( A, b, s, M, tol, maxit );

%
% RITZ_IDRS: Solve linear system with IDRS and compute s ritzvalues and ritzvectors 
% using the IDR Hessenberg matrix. Using this spectral information an initial search 
% space and/or parameters omega are computed for subsequent solves
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen
%

% Inital IDRS solve
   [x, flag, relres, iter, resvec, H] = idrs( A, b, s, M, tol, maxit );
   n = length(b);
   [m1,m] = size(H);
%
% Compute ritzvalues
   [ritzval ritzvec_s] = select( H(1:m,1:m), 'SM' );
   ritzvec_s = ritzvec_s(:,1:s);

   U0 = zeros(n,s);
%
% Compute the eigenvectors
   W = zeros(n,s+1);
   w = b;
   ind = zeros(m,1);
   k = 0;
   for i = 1:m
      k = k + 1;
      if ( k > s+1 ) k = 1; end;
      ind(i) = k;
      W(:,k) = w;
      U0 = U0 + w*ritzvec_s(i,:);
      Aw = A*w;
      low = max(1,i-s);
      w = (Aw - W(:,ind(low:i))*H(low:i,i)) /H(i+1,i);
   end
%
   omega = 1./leja(ritzval);

return
