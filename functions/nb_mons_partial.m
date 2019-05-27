function nb=nb_mons_partial(n,d),
%NB_MONS_PARTIAL   Return the number of monomials in a partial monomial basis.
%  
% SIGNATURE
% nb=nb_mons_partial(nvar,d);
%
% DESCRIPTION
% Returns the number of monomials in a partial monomial basis of nvar unknowns 
% and total degree d (partial: only degree d) 
%
% INPUTS
%    nvar   =   number of variables in input system
%    d      =   degree of partial monomial basis
%
% OUTPUTS
%    nb     =   number of monomials in partial monomial basis of nvar variables and degree d
%
% EXAMPLE
%
% CALLS
%    
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

nb = (factorial(d+n-1)/factorial(d))/factorial(n-1);

end
