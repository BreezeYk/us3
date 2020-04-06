function [P,l] = pnseq (pu,varargin)

% PNSEQ PN codes.
%    [P,l] = PNSEQ(PU) generates a pseudonoise (PN) sequence P from a
%    binary linear feedback shift register.  l gives the length of the
%    code.
%    The shift register is described by the generator polynomial PU.
%
%    PNSEQ(...,RU) generates a PN sequences given RU, the initial state of
%    the shift register.  RU must be a vector of the same length as the
%    generator polynomial.
%
%    Note that a {0,1}-valued binary sequence corresponds to a
%    cos(pi*{0,1}) = {1,-1}-transmitted sequence.

% Joe Henning - Fall 2006

P = 0;
l = 0;

error(nargchk(1,2,nargin));
[ru,n,msg] = parseinput(pu,varargin{:});
error(msg);

l = 2^n - 1;

P = zeros(1,l);
for k = 1:l
   temp = mod( fliplr(pu)*ru.', 2 );
   P(k) = ru(n);
   ru = [temp ru(1:n-1)];
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
