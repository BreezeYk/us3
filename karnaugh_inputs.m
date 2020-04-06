function K = karnaugh_inputs (n)

% KARNAUGH_INPUTS Boolean inputs for a Karnaugh map.
%    K = KARNAUGH_INPUTS(N) produces a table (array) of boolean inputs for
%    a Karnaugh map.  Each boolean input has N bits, and there are 2^N
%    entries in the table.

% Joe Henning - Fall 2006

A = [0; 1];
K = zeros(2^n,n);
K(:,n) = repmat(A,2^(n-1),1);
for k = n-1:-1:2
   K(:,k) = repmat([repmat(A(1),2^(n-k),1); repmat(A(2),2^(n-k),1)], 2^(n-(n-k+1)), 1);
end
K(:,1) = [repmat(A(1),2^(n-1),1); repmat(A(2),2^(n-1),1)];
