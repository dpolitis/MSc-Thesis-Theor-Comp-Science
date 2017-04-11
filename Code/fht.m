function fOut=fht(fIn)

%**************************************************************************
%
%  fht computes a radix-2 1D Fast Hartley Transform.
%
%  Discussion:
%
%    The discrete Hartley transform h of a set of data a is
%
%      h(i) = 1/sqrt(N) * sum (0<=j<=N-1) a(j) * cas(2*pi*i*j/N)
%
%    The data and coefficients are indexed from 1 to N.
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
%  Remarks:
%
%    Length of fIn() Should be 2^p. If it is not, the code pads input
%    with zeros and removes them at end of transform.
%
%    Sequential FHT Algorithm
%    Input is bit-reversed in Code in fIn(1..N)
%    Output in Normal Order in f(1..N)
%
%**************************************************************************

[k,l]=size(fIn); % Check if transpose is needed
if l>k
    fIn=fIn'; n=l; flag=1;
else
    n=k; flag=0;
end;

[f,ldn]=log2(n);
Nold=0; % To be used later as flag

if (f~=0.5)  % Pad array with zeros so that length is power of 2
    fIn=[fIn;zeros(2^ldn-n,1)];
    Nold=n;  % Keep to cut off redundant values at end of transform
    n=2^ldn;
end;

fIn=bitrevorder(fIn); % Reverse bit order

for ldm=1:ldn
    m=2^ldm;
    mh=m/2;
    m4=m/4;
    arg=pi/mh;
    
    for r=0:m:n-m
        for j=1:m4-1 % Hartley shift
            k=mh-j;
            
            u=fIn(r+mh+j+1);
            v=fIn(r+mh+k+1);
            
            c=cos(j*arg);
            s=sin(j*arg);

            fIn(r+mh+j+1)=c*u+s*v;
            fIn(r+mh+k+1)=s*u-c*v;
        end
        
        for j=0:mh-1
            u=fIn(r+j+1);
            v=fIn(r+j+mh+1);
            
            fIn(r+j+1)=u+v;
            fIn(r+j+mh+1)=u-v;
        end
    end
end

% Truncate last values if original data length was not power of 2
if (Nold~=0)
    fIn=fIn(1:Nold);
    n=Nold;
end

if flag==1; % Transpose if needed
    fIn=fIn';
end;

fOut=fIn/sqrt(n); % Make DHT it's own inverse.
end
