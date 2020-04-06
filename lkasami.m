function [Kl,l] = lkasami (pu,pv,varargin)

% LKASAMI Large-set of Kasami codes.
%    [Kl,l] = LKASAMI(PU,PV) generates the large set of Kasami sequences Kl
%    from two binary linear feedback shift registers.  l gives the length
%    of the code.
%    The shift registers are described by the preferred polynomials PU and
%    PV.  PU and PV must be vectors of the same length.
%
%    LKASAMI(...,RU,RV) generates the large set of Kasami sequences given
%    RU and RV, the initial states of the U and V shift registers.  RU and
%    RV must be vectors of the same length as the preferred polynomials.
%
%    Note that a {0,1}-valued binary sequence corresponds to a
%    cos(pi*{0,1}) = {1,-1}-transmitted sequence.

% Joe Henning - Fall 2006

Kl = 0;
l = 0;

error(nargchk(2,4,nargin));
[ru,rv,n,msg] = parseinput(pu,pv,varargin{:});
error(msg);

l = 2^n - 1;

t = 2^(n/2);
t2 = 2^((n+2)/2);

%fprintf('Cross-correlation can take the values: %.3f, %.3f, %.3f., %.3f, %.3f\n',-1,-t-1,t-1,-t2-1,t2-1);

[g,G,l,OG] = gold(pu,pv,0,ru,rv);

u = G(1,:);
v = G(2,:);

q = 2^(n/2) + 1;

w = downsample(u,q);

Kl = G;
for k = 1:size(G,1)
   for m = 0:2^(n/2)-2
      shifted_w = circshift(w,[0 -m]);
      upw = repmat( shifted_w, 1, q );
      Kl = [Kl; mod( G(k,:) + upw, 2 )];
   end
end



%--------------------------------------------------------------------------
function [ru,rv,n,msg] = parseinput(pu,pv,varargin)
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

if ( mod(n,2) ~= 0 )
   msg = 'n must be even.';
   return
end

% Assume no register values until proven otherwise
found_ru = 0;
found_rv = 0;
ru = ones(1,n);
rv = ones(1,n);

switch nargin
   case 3
      % is ru
      if isnumeric(varargin{1})
         if length(varargin{1}) == n
            found_ru = 1;
            ru = varargin{1};
         end
      end
   case 4
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
