function [C] = corr (X,varargin)
% CORR Periodic and aperiodic cross-correlation function estimates.
%    C = CORR(A,B), where A and B are length M vectors (M>1), returns
%    the length 2*M-1 aperiodic cross-correlation sequence C.  If A
%    and B are of different length, the shortest one is zero-padded.
%
%    CORR(A), when A is a vector, is the aperiodic autocorrelation
%    sequence.
%
%    CORR(..,CORRTYPE) specifies the type of cross-correlation according
%    to CORRTYPE:
%       'a' - aperiodic cross-correlation (this is the default)
%       'p' - periodic cross-correlation
%
%    A periodic cross-correlation sequence C is of length M.
% Joe Henning - Fall 2006
C = 0;
error(nargchk(1,3,nargin));
[autoFlag,Y,corrType,msg] = parseinput(X,varargin{:});
error(msg);
Nx = length(X);
Ny = length(Y);
if Nx < Ny
   X = [X zeros(1,Ny-Nx)];
   N = Ny;
elseif Ny < Nx
   Y = [Y zeros(1,Nx-Ny)];
   N = Nx;
else
   N = Nx;
end
if ( strcmp(corrType,'a') | strcmp(corrType,'A') )
   %C = xcorr(X);
   C = zeros(1,2*N-1);
   for k = 1:2*N-1
      for l = 1:1:N
         if ( l-k+N > 0 & l-k+N <= N )
            C(k) = C(k) + X(l)*conj(Y(l-k+N));
         end
      end
   end
else   % corrType = 'p' | corrType = 'P'
   C = zeros(1,N);
   for l = 0:1:N-1
      for n = 0:1:N-1
         yii = n+l+1;
         if ( yii > N )
            yii = mod( yii, N );
         end
         C(l+1) = C(l+1) + X(n+1)*conj(Y(yii));
      end
   end
end
%--------------------------------------------------------------------------
function [autoFlag,Y,corrType,msg] = parseinput(X,varargin)
%   Parse the input and determine optional parameters:
%
%   Outputs:
%      autoFlag  - 1 if autocorrelation, 0 if xcorrelation
%      Y         - vector Y
%      corrType  - string with the type of correlation wanted
%      msg       - possible error message
% Set some defaults:
% Assume aperiodic autocorrelation until proven otherwise
msg = '';
corrType = 'a';
autoFlag = 1;
Y = X;
errMsg = 'Input argument is not recognized.';
switch nargin
   case 2
      % Can be (x,y) or (x,corrType)
      if ischar(varargin{1})
         % Second arg is corrType
         corrType = varargin{1};
      elseif isnumeric(varargin{1})
         autoFlag = 0;
         Y = varargin{1};
      else
         % Not recognized
         msg = errMsg;
         return
      end
   case 3
      % Must be (x,y,corrType)
      autoFlag = 0;
      Y = varargin{1};
      corrType = varargin{2};
end
[mx,nx] = size(X);
if ( (mx ~= 1 & nx ~= 1) | (mx == 1 & nx == 1) )
   msg = 'A must be a vector.';
   return
end
if ~autoFlag
   [my,ny] = size(Y);
   if ( (my ~= 1 & ny ~= 1) | (my == 1 & ny == 1) )
      msg = 'B must be a vector.';
      return
   end
end
if ( ~strcmp(corrType,'a') & ~strcmp(corrType,'A') & ~strcmp(corrType,'p') & ~strcmp(corrType,'P') )
   fprintf('Unknown corrType.  corrType set to ''a''.\n');
   corrType = 'a';
end