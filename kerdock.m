function [k,v,l] = kerdock (pu,varargin)

% KERDOCK Kerdock codes.
%    [k,v,l] = KERDOCK(PU) generates a binary Kerdock code k and
%    its corresponding meander sequence v.  l gives the length of the
%    code.
%    The shift register is described by a quaternary characteristic
%    polynomial PU.
%
%    KERDOCK(...,RU) generates a Kerdock code given RU, the initial
%    state of the shift register.  Elements of RU should be in Z4, and
%    RU must be a vector whose length is one less than the size of PU
%    (i.e. a trailing binary 1 is not assumed with the characteristic
%    polynomial).  Changing the register's initial loading results in
%    2^(n-1) different Kerdock codes.
%
%    Note that a {0,1}-valued binary sequence corresponds to a
%    cos(pi*{0,1}) = {1,-1}-transmitted sequence.
%
%    The following table provides a short list of quaternary
%    characteristic polynomials:
%       n    Characteristic Polynomial
%       --   -------------------------
%       5    [ 1 0 0 1 2 1]
%       7    [ 1 0 0 2 0 0 1 1]
%       9    [ 1 0 0 0 0 1 0 2 0 1]
%       11   [ 1 0 0 0 0 0 0 0 0 1 2 1]
%       13   [ 1 3 1 0 0 0 2 0 0 0 0 0 1 1]
%       15   [ 1 0 0 0 0 0 0 2 0 0 0 0 0 0 1 1]
%       17   [ 1 0 0 0 0 0 0 2 0 0 0 0 0 0 1 0 0 1]

% Joe Henning - Jan 2012

k = 0;
v = 0;
l = 0;

error(nargchk(1,2,nargin));
[ru,n,msg] = parseinput(pu,varargin{:});
error(msg);

l = 2^n - 1;

qu = mod(-pu(2:end), 4);

k = zeros(1,l);
for m = 1:2*l
   temp = mod( qu*ru.', 4 );
   k(m) = ru(n);
   ru = [temp ru(1:n-1)];
end

% MSB mapping
i = find(k <= 1);
k(i) = 0;
i = find(k >= 2);
k(i) = 1;

% Generate corresponding sequence through bit meander
n = 0:length(k)-1;
v = mod(mod(n,2) + k, 2);



%--------------------------------------------------------------------------
function [ru,n,msg] = parseinput(pu,varargin)
%    Parse the input and determine optional parameters:
%
%    Outputs:
%       ru  - initial values of register u
%       n   - number of register bits
%       msg - possible error message

% Set some defaults:
ru = 0;
n = 0;
msg = '';

[mu,nu] = size(pu);
if ( (mu ~= 1 & nu ~= 1) | (mu == 1 & nu == 1) )
   msg = 'pu must be a vector.';
   return
end

n = length(pu) - 1;

if ( mod(n,2) ~= 1 )
   msg = 'n must be odd.';
   return
end

% Assume no register values until proven otherwise
found_ru = 0;
ru = ones(1,n);

switch nargin
   case 2
      % is ru
      if isnumeric(varargin{1})
         if length(varargin{1}) == n
            found_ru = 1;
            ru = varargin{1};
         end
      end
end

if ~found_ru
   fprintf('Using default ru.\n');
end
