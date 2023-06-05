%
% Internal subroutine for the idrs package
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen

function om = comp_om( t, r, kappa )

   if ( kappa == 0 )
      om = dot(t,r)/dot(t,t);
   else
      nr = norm(r);
      nt = norm(t);
      tr = dot(t,r);
      rho = abs(tr/(nt*nr));
      om=tr/(nt*nt);
      if ( rho < kappa )
         om = om*kappa/rho;
      end
   end

return
