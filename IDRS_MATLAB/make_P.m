%
% Internal subroutine for the idrs package
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen

function P = make_P( r, s, method )

n = length(r);
if ( method == 1 )
   rng('default');
   P = randn(n,s);
   P = orth(P);
else
   P = r;
end 

return
