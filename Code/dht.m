function fOut=dht(fIn)

%**************************************************************************
%
%  dht computes a "slow" 1-D Hartley Transform.
%
%  Discussion:
%
%    The discrete Hartley transform h of a set of data a is
%
%      h(i) = 1/sqrt(N) * sum (0<=j<=N-1) a(j) * cas(2*pi*i*j/N)
%
%    The data and coefficients are indexed from 0 to N-1.
%
%    With the above normalization factor of 1/sqrt(N), the Hartley
%    transform is its own inverse.
%
%    This routine is provided for illustration and testing.  It is
%    inefficient relative to optimized routines.
%
%  Licensing:
%
%    This code is distributed under the GNU GPLv3 license.
%    Copy of the GPLv3 License can be found in the following URL:
%    http://www.gnu.org/licenses/gpl-3.0.html
%
%  Modified:
%
%    10 Nov 2011
%
%  Author:
%
%    Dimitrios Politis
%
%  Reference:
%
%    Ralph Hartley,
%    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
%    Proceedings of the Institute of Radio Engineers,
%    Volume 30, pages 144-150, 1942.
%
%  Parameters:
%
%    Input, real vector fIn with rows, the data to be transformed.
%
%    Output, real vector fOut, the transformed data.
%
%**************************************************************************

[m,n]=size(fIn); % Check if transpose is needed
if n>m
    fIn=fIn'; N=n; flag=1;
else
    N=m; flag=0;
end;

fOut=zeros(N,1);
arg0=2*pi/N;
arg1=pi/4;
arg2=sqrt(2);
n=(0:N-1)';

for k=0:N-1 % Begin 1-D Hartley Transform
   fOut(k+1)=sum(fIn.*(sin(arg0*n*k+arg1))); % cas(x)=sqrt(2)*sin(x+pi/4)
end         % End 1-D Hartley Transform

if flag==1;
    fOut=fOut';
end;

fOut=arg2*fOut./sqrt(N);
end