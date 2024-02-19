function [eival, eivec] = select( H, target );
%
% Compute eigenvalues and eigenvector of small problem
% Select the ones closest to the target
%

eival = [];
eivec = [];
if ( isempty(H) ) return; end;

comp_eivec = ( nargout == 2 );

if comp_eivec
   [V,D] = eig(H);
   ei = diag(D);
else
   ei = eig(H);
end
n = length(ei);

if ( strcmp(target,'SM') || strcmp(target,'sm') ) 
   [dum, ind] = sort(abs(ei));
elseif ( strcmp(target,'LM') || strcmp(target,'lm') )
   [dum, ind] = sort(-abs(ei));
elseif ( strcmp(target,'SR') || strcmp(target,'sr') )
   [dum, ind] = sort(real(ei));
elseif ( strcmp(target,'LR') || strcmp(target,'lr') )
   [dum, ind] = sort(-real(ei));
elseif ( strcmp(target,'SI') || strcmp(target,'si') )
   [dum, ind] = sort(imag(ei));
elseif ( strcmp(target,'LI') || strcmp(target,'li') )
   [dum, ind] = sort(-imag(ei));
else
   error('Illegal value for TARGET');
end

eival = ei(ind);
if comp_eivec
   eivec = V(:,ind);
end

return
