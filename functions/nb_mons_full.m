function nb = nb_mons_full(n,d),
%NB_MONS_FULL   Return the number of monomials in a full monomial basis.
%  
% SIGNATURE
% nb=nb_mons_full(nvar,d);
%
% DESCRIPTION
% Returns the number of monomials in a full monomial basis of nvar unknowns 
% and total degree d (full: degrees 0:d)
%
% INPUTS
%    nvar   =   number of variables in input system
%    d      =   maximal total degree of full monomial basis
%
% OUTPUTS
%    nb     =   number of monomials in full monomial basis of nvar variables and degree d
%
% EXAMPLE
%
% CALLS
%   nb_mons_partial 
%
% AUTHOR
%   Philippe Dreesen (philippe.dreesen@gmail.com)
%   June 2010
%

nb=0;
for i=0:d,
	nb = nb + nb_mons_partial(n,i);
end
