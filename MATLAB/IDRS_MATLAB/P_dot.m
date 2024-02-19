%
% Internal subroutine for the idrs package
%
% This software is distributed under the MIT License:
% http://www.opensource.org/licenses/mit-license.php
% Copyright:(c) 2023 Martin van Gijzen

function f = P_dot( P, g, s )

f = zeros(s,1);
for i = 1:s
    f(i) = dot(P(:,i),g);
end

return
