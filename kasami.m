function [Ks,l] = kasami (pu,varargin)

% KASAMI Small-set of Kasami codes.
%    [Ks,l] = KASAMI(PU) generates the small set of Kasami sequences Ks
%    from a binary linear feedback shift register.  l gives the length of
%    the code.
%    The shift register is described by the generator polynomial PU.
%
%    KASAMI(...,RU) generates the small set of Kasami sequences given RU,
%    the initial state of the U shift register.  RU must be a vector of the
%    same length as the generator polynomial.
%
%    Note that a {0,1}-valued binary sequence corresponds to a
%    cos(pi*{0,1}) = {1,-1}-transmitted sequence.

% Joe Henning - Fall 2006

Ks = 0;
l = 0;

error(nargchk(1,2,nargin));
[ru,n,msg] = parseinput(pu,varargin{:});
error(msg);

l = 2^n - 1;

t = 2^(n/2);

%fprintf('Cross-correlation can take the values: %.3f, %.3f, %.3f.\n',-1,-t-1,t-1);

u = zeros(1,l);
for k = 1:l
   temp = mod( fliplr(pu)*ru.', 2 );
   u(k) = ru(n);
   ru = [temp ru(1:n-1)];
end

q = 2^(n/2) + 1;

w = downsample(u,q);

Ks = u;
%Ks = [];
for k = 0:2^(n/2)-2
   shifted_w = circshift(w,[0 -k]);
   upw = repmat( shifted_w, 1, q );
   Ks = [Ks;  mod( u + upw, 2 )];
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
