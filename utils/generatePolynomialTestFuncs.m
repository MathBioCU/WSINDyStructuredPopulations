function [testFuncs, testDerivs,sups] = generatePolynomialTestFuncs(a, b, p, q)
% generatePolynomialTestFuncs - Generates test functions of the form
% f(t) = C(t - a)^p (b - t)^q when a < t < b and 0 otherwise, along with
% their derivatives.
% Inputs:
%   a, b       - Vectors of endpoints of the supports of the test function ( a_i < t < b_i).
%   p          - Power for (t - a).
%   q          - Power for (b - t).
% Outputs:
%   testFuncs  - Cell array of function handles for the test functions.
%   testDerivs - Cell array of function handles for derivatives of test 
%                functions.
%   sups       - Cell array of vectors specifiying the supports of test
%                functions given by inputs a and b.



% Initialize cell arrays to store test functions and their derivatives
testFuncs = cell(length(a), 1);
testDerivs = cell(length(a), 1);
sups = cell(length(a),1);


% Generate test functions and their derivatives for each combination of p and q
for i =1:length(a)
               
        C = 1/(p^p *q^q) * ((p+q)/(b(i)-a(i)))^(p+q);  % Normalization constant
        testFuncs{i} = @(t)  C * (t > a(i) & t < b(i)) .* (t - a(i)).^p .* (b(i) - t).^q;   
        testDerivs{i}= @(t) C * (t > a(i) & t < b(i)) .* ( p*(t - a(i)).^(p-1) .* (b(i) - t).^q-q*(t - a(i)).^(p) .* (b(i) - t).^(q-1));

               
        sups{i}=[a(i),b(i)];
end

end
