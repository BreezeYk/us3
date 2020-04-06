function [g,G,l,OG] = gold (pu,pv,shift,varargin)

% GOLD Gold codes.
%    [g,G,l,OG] = GOLD(PU,PV,SHIFT) generates a binary Gold code g and Gold
%    code sequence G from two binary linear feedback shift registers.  l
%    gives the length of the code.  OG is an orthogonal Gold code
%    sequence.
%    The shift registers are described by the preferred polynomials PU and
%    PV.  PU and PV must be vectors of the same length.  SHIFT is an
%    integer scalar which specifies the offset of the singular Gold code
%    from the initial time.
%
%    GOLD(..,RU,RV) generates Gold codes given RU and RV, the initial
%    states of the U and V shift registers.  RU and RV must be vectors of
%    the same length as the preferred polynomials.
%
%    Note that a {0,1}-valued binary sequence corresponds to a
%    cos(pi*{0,1}) = {1,-1}-transmitted sequence.
%
%    The following table provides a short list of preferred pairs
%    (a trailing binary 1 is always assumed):
%       n    N      Preferred Polynomial 1   Preferred Polynomial 2
%       --   ----   ----------------------   ----------------------
%       5    31     [ 5 2 ] = [1 0 0 1 0]    [ 5 4 3 2 ] = [1 1 1 1 0]
%       6    63     [ 6 1 ] = [1 0 0 0 0 1]  [ 6 5 2 1 ] = [1 1 0 0 1 1]
%       7    127    [ 7 3 ] = ...            [ 7 3 2 1 ] = ...
%       9    511    [ 9 4 ] = ...            [ 9 6 4 3 ] = ...
%       10   1023   [ 10 3 ] = ...           [ 10 8 3 2 ] = ...
%       11   2047   [ 11 2 ] = ...           [ 11 8 5 2 ] = ...
%    See http://www.cs.umbc.edu/~lomonaco/f97/442/Peterson_Table.html for a
%    table of irreducible polynomials over GF(2).
%
%    Given a polynomial representation in octal, the polynomial in binary
%    is de2bi(oct2dec(x),'left-msb').
%    For example, octal [23] is [1 0 0 1 1] in binary.  Since a trailing
%    binary 1 is always assumed, the preferred polynomial is [1 0 0 1].

% Joe Henning - Fall 2006

g = 0;
G = 0;
l = 0;
O = 0;

error(nargchk(3,5,nargin));
[ru,rv,n,msg] = parseinput(pu,pv,shift,varargin{:});
error(msg);

l = 2^n - 1;

if mod(n,2) == 1   % odd
   t = 1 + 2^((n+1)/2);
else   % even
   t = 1 + 2^((n+2)/2);
end

%fprintf('Cross-correlation can take the values: %.3f, %.3f, %.3f.\n',-t,-1,t-2);

u = zeros(1,l);
for k = 1:l
   temp = mod( fliplr(pu)*ru.', 2 );
   u(k) = ru(n);
   ru = [temp ru(1:n-1)];
end

v = zeros(1,l);
for k = 1:l
   temp = mod( fliplr(pv)*rv.', 2 );
   v(k) = rv(n);
   rv = [temp rv(1:n-1)];
end

shifted_v = circshift(v,[0 -shift]);

g = mod( u + shifted_v, 2 );

G = [u; v];
%G = [];
for k = 0:l-1
   shifted_v = circshift(v,[0 -k]);
   G = [G; mod( u + shifted_v, 2 )];
end

OG = [G zeros(size(G,1),1)];



%--------------------------------------------------------------------------
function [ru,rv,n,msg] = parseinput(pu,pv,shift,varargin)
%    Parse the input and determine optional parameters:
%
%    Outputs:
%       ru  - initial values of register u
%       rv  - initial values of register v
%       n   - number of register bits
%       msg - possible error message

% Set some defaults:
ru = 0;
rv = 0;
n = 0;
msg = '';

[mu,nu] = size(pu);
if ( (mu ~= 1 & nu ~= 1) | (mu == 1 & nu == 1) )
   msg = 'pu must be a vector.';
   return
end

[mv,nv] = size(pv);
if ( (mv ~= 1 & nv ~= 1) | (mv == 1 & nv == 1) )
   msg = 'pv must be a vector.';
   return;
end

if ( length(pu) ~= length(pv) )
   msg = 'pu and pv must be of the same length.';
   return
end

n = length(pu);

if ~isnumeric(shift)
   msg = 'shift must be an integer.';
   return
end

% Assume no register values until proven otherwise
found_ru = 0;
found_rv = 0;
ru = ones(1,n);
rv = ones(1,n);

switch nargin
   case 4
      % is ru
      if isnumeric(varargin{1})
         if length(varargin{1}) == n
            found_ru = 1;
            ru = varargin{1};
         end
      end
   case 5
      % is ru and rv
      if isnumeric(varargin{1})
         if length(varargin{1}) == n
            found_ru = 1;
            ru = varargin{1};
         end
      end
      if isnumeric(varargin{2})
         if length(varargin{2}) == n
            found_rv = 1;
            rv = varargin{2};
         end
      end
end

if ~found_ru
   fprintf('Using default ru.\n');
end

if ~found_rv
   fprintf('Using default rv.\n');
end
