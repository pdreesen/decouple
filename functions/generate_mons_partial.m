function base = generate_mons_partial(n,d)
%GENERATE_MONS_PARTIAL Generate basis of monomials of given degree.
%  
% SIGNATURE
% base = generate_mons_partial(nvar,d)
%
% DESCRIPTION
% Generate a partial basis (degree 'd'; 'n' variables) of monomials as a 
% matrix. Each row of 'base' refers to a n-tuple of exponents of a monomial 
% where each column corresponds to a variable.
% 
% INPUTS
%    nvar       =    number of variables
%    d          =    desired degree
%
% OUTPUTS
%    base       =    partial basis of monomials 
%
% EXAMPLE
%
% CALLS
%    [none]
%
% AUTHOR
%    Philippe Dreesen (philippe.dreesen@gmail.com)
%    June 2010
%

% TODO: preallocate base and work with index of writeatrow or something

if n <= 1
    base = d;
else
    base = [];
    for i = d :-1: 0
        temp = generate_mons_partial(n-1,d-i);
        base = [base; i*ones(size(temp,1),1) generate_mons_partial(n-1,d-i)];
    end
end
end
