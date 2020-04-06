function [w,W] = weil (p,varargin)

% WEIL Weil codes.
%    [w,W] = WEIL(P) generates a binary Weil code w and Weil code
%    sequence W.  P specifies the length of the code; it must be
%    a prime number.
%
%    WEIL(...,CODENUM) returns the NUM-th Weil code.  If missing,
%    default is CODENUM = 1.

% Joe Henning - 5 Jan 2012

w = [];
W = [];

error(nargchk(1,2,nargin));
[num,msg] = parseinput(p,varargin{:});
error(msg);

t = 1 + 2*sqrt(p);

%fprintf('Cross-correlation has a maximum value of %.3f\n',t);

for k = 1:(p-1)/2
   count = 1;
   for i = 0:p-1
      if (i+k > p-1)
         W(k,count) = mod(leg(i,p) + leg(i+k-p,p),2);
      else
         W(k,count) = mod(leg(i,p) + leg(i+k,p),2);
      end
      count = count + 1;
   end
end

w = W(num,:);



%--------------------------------------------------------------------------
function [num,msg] = parseinput(p,varargin);

% Outputs: p   - code length
%          num - desired code index
%          msg - error message

num = 1;   % Assume first code
msg = [];

if (p < 3)
   msg = sprintf('Error:  p must be greater than 2.');
   return
end

if ~isprime(p)
   msg = sprintf('Error:  p must be prime.');
   return
end

if nargin > 1
   if isnumeric(varargin{1})
      num = varargin{1};
   else
      msg = sprintf('Error:  num must be an integer scalar.');
      return
   end
end

num = fix(num);

if (num < 1) | (num > (p-1)/2)
   msg = sprintf('Error:  num must be >= 1 and <= (p-1)/2.');
   return
end



%--------------------------------------------------------------------------
function L = leg(i,p)

% Compute the Legendre sequence

if (i == 0)
   L = 1;
else
   if (mod(i,p) == 0)
      L = NaN;   % p divides i (doesn't happen for primes)
   else
      squaremodp = 0;
      z = 1;
      while (z <= (p-1)/2)
         if (mod(z*z,p) == i)
            squaremodp = 1;
            break;
         end
         z = z + 1;
      end

      if (squaremodp)
         L = 0;   % i is a square(mod p)
      else
         L = 1;   % i is not a square(mod p)
      end
   end
end