function fOut = cfht2fft(fIn,fft_s)

%**************************************************************************
%
%  cfht2fft computes a complex 1D-FFT or 1D-IFFT through
%  1-D Hartley Transform. Forward or Inverse transform depends on fft_s.
%  Use -1 for forward and 1 for inverse transform.
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
%    Input, complex vector fIn with rows, the data to be transformed.
%
%    Input, Integer fft_s, the Fourrier transform sign
%
%    Output, complex vector fOut, the transformed data.
%
%**************************************************************************

[m,n]=size(fIn); % Check if transpose is needed
if n>m
    fIn=fIn'; flag=1; N=n;
else
    flag=0; N=m;
end;

% Compute FHT for real and imaginary parts
temp(:,1)=fht(real(fIn));
temp(:,3)=fht(imag(fIn));
temp(:,2)=[temp(1,1);flipud(temp(2:end,1))];
temp(:,4)=[temp(1,3);flipud(temp(2:end,3))];

realPart=temp(:,1)+temp(:,2)-fft_s*(temp(:,3)-temp(:,4));
imPart=temp(:,3)+temp(:,4)+fft_s*(temp(:,1)-temp(:,2));

% Comply with Matlab FFT function
fOut=sqrt(N)*0.5*(realPart+1i*imPart);

if flag==1; % Transpose
    fOut=fOut';
end;

end