function [Hq,l,q_set] = goldlike (pu,varargin)

% GOLDLIKE Gold-like and Dual-BCH codes.
%    [Hq,l,q] = GOLDLIKE(PU) generates a Gold-like or Dual-BCH sequence
%    Hq from a binary linear feedback shift register.  l gives the length
%    of the code.  q_set is a set of integers such that gcd(q,2^n-1) = 3.
%    The first member of q_set is used in the algorithm.
%
%    GOLDLIKE(...,RU) generates Gold-like or Dual-BCH sequence given RU,
%    the initial state of the U shift register.  RU must be a vector of the
%    same length as the generator polynomial.
%
%    Note that a {0,1}-valued binary sequence corresponds to a
%    cos(pi*{0,1}) = {1,-1}-transmitted sequence.

% Joe Henning - Fall 2006

Hq = 0;
l = 0;

error(nargchk(1,2,nargin));
[ru,n,msg] = parseinput(pu,varargin{:});
error(msg);

l = 2^n - 1;

u = zeros(1,l);
for k = 1:l
   temp = mod( fliplr(pu)*ru.', 2 );
   u(k) = ru(n);
   ru = [temp ru(1:n-1)];
end

q_set = [];
for k = 1:l
   if ( gcd(k,l) == 3 )% & ( mod(l,k) == 0 )
      q_set = [q_set k];
   end
end

q = q_set(1);

lprime = l/q;

Hq = u;
%Hq = [];
for k = 0:1:2
   shifted_u = circshift(u,[0 -k]);
   v = downsample( shifted_u, q );
   for m = 0:lprime-1
      shifted_v = circshift(v,[0 -m]);
      upv = repmat( shifted_v, 1, q );
      Hq = [Hq; mod( u + upv, 2 )];
   end
end



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

n = length(pu);

if ( mod(n,2) ~= 0 )
   msg = 'n must be even.';
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
