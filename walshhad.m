function [w,W,h,H] = walshhad (l,varargin)

% WALSHHAD Walsh and Hadamard codes.
%    [w,W,h,H] = WALSHHAD(L) generates a binary Walsh code w, Walsh
%    code sequence W, binary Hadamard code h, and Hadamard code sequence
%    H.  L specifies the length of the code.
%
%    WALSHHAD(...,CODENUM) returns the NUM-th Walsh/Hadamard code.  If
%    missing, default is CODENUM = 1.

% Joe Henning - 13 Jan 2008

w = [];
W = [];
h = [];
H = [];

error(nargchk(1,2,nargin));
[num,msg] = parseinput(l,varargin{:});
error(msg);

H = hadamard(l);
h = H(num,:);

for k = 1:length(h)
   zc = 1;   % matlab starts array indices at 1
   for m = 2:length(h)
      if (H(k,m) ~= H(k,m-1))
         zc = zc + 1;
      end
   end
   W(zc,:) = H(k,:);
end

w = W(num,:);



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
