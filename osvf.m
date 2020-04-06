function [o,O] = osvf (l,varargin)

% OSVF Orthogonal variable spreading factor (OSVF) codes.
%    [o,O] = OSVF(L) generates a binary OSVF code and OSVF code sequence
%    O.  L specifies the length of the code (and spreading factor).
%
%    OSVF(...,CODENUM) returns the NUM-th OSVF code.  If missing, default
%    is CODENUM = 1.

% Joe Henning - 13 Jan 2008

o = [];
O = [];

error(nargchk(1,2,nargin));
[num,msg] = parseinput(l,varargin{:});
error(msg);

C_old = 1;
if (l > 1)
   n = 1;
   while 2*n <= l
      n = 2*n;
      C_new = [];
      for k = 1:n/2
         C_new = [C_new; C_old(k,:) C_old(k,:); C_old(k,:) -C_old(k,:)];
      end
      C_old = C_new;
   end
end

o = C_old(num,:);
O = C_old;



%--------------------------------------------------------------------------
function [num,msg] = parseinput(l,varargin);

% Outputs: l   - code length
%          num - desired code index
%          msg - error message

N = [];
num = 1;   % Assume first code
msg = [];

if (l < 1)
   msg = sprintf('Error:  l must be an integer scalar greater than 0.');
   return
end

N = log2(l);

if (N - fix(N) > 1E-8)
   msg = sprintf('Error:  l must be a power of 2.');
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

if (num < 1) | (num > l)
   msg = sprintf('Error:  num must be >= 1 and <= l.');
   return
end
