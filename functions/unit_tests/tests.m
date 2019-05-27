% test function for 'constructVariablesPolynomial'

base = [2; 3];
result = constructVariablesPolynomial(base, 5, 0);
resultAnswer = [1 2 4 8 16 32; 1 3 9 27 81 243];
assert(isequal(result,resultAnswer))

base = [2; 3];
result = constructVariablesPolynomial(base, 5, 1);
resultAnswer = [0 1 4 12 32 80; 0 1 6 27 108 405];
assert(isequal(result,resultAnswer))

base = [2; 3];
result = constructVariablesPolynomial(base, 5, 2);
resultAnswer = [0 0 2 12 48 160; 0 0 2 18 108 540];
assert(isequal(result,resultAnswer))

base = sym('u', [1 r]).';
result = constructVariablesPolynomial(base, 5, 2);
resultAnswer = [ 0, 0, 2, 6*base(1), 12*base(1)^2, 20*base(1)^3; 0, 0, 2, 6*base(2), 12*base(2)^2, 20*base(2)^3];
assert(isequal(result,resultAnswer))