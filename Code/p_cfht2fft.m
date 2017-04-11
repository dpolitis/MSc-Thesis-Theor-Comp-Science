function fOut = p_cfht2fft(fIn,fft_s)

%**************************************************************************
%
%  p_cfht2fft computes a complex 1D-FFT or 1D-IFFT through parallel
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
if (labindex==1) % Lab1 is master node
    [m,n]=size(fIn); % Check if transpose is needed
    if n>m
        fIn=fIn'; flag=1; N=n;
    else
        flag=0; N=m;
    end;

    labSend(fIn,2,0); % Send input vector to Lab 2
    labSend(fft_s,2,1); % Send fourrier transform sign to Lab2
    
    re(:,1)=fht(real(fIn)); % Compute FHT for real part
    re(:,2)=[re(1,1);flipud(re(2:end,1))]; % fht(N-k)
    
    re(:,3)=re(:,1)+re(:,2); re(:,4)=re(:,1)-re(:,2);
    
    labSend(re(:,3:4),2,2); % Send real sumdiffs to Lab2
    im=labReceive(2,3); % Receive imaginary sumdiffs from Lab2
    
    realPart=re(:,3)-fft_s*im(:,2);

    imPart=labReceive(2,4); % Receive imaginary part from Lab2
        
    % Comply with Matlab FFT function
    fOut=sqrt(N)*0.5*(realPart+1i*imPart);

    if flag==1; % Transpose
        fOut=fOut';
    end;
end;

if (labindex==2) % Lab2 executes from now on
    fIn=labReceive(1,0);
    fft_s=labReceive(1,1); % Receive fourrier transform sign from Lab1
    
    im(:,1)=fht(imag(fIn)); % Compute FHT for imaginary part
    im(:,2)=[im(1,1);flipud(im(2:end,1))]; % fht(N-k)

    im(:,3)=im(:,1)+im(:,2); im(:,4)=im(:,1)-im(:,2);
    
    re=labReceive(1,2); % Receive real sumdiffs from Lab1
    labSend(im(:,3:4),1,3); % Send imaginary sumdiffs to Lab1
    
    imPart=im(:,1)+fft_s*re(:,2);

    labSend(imPart,1,4); % Send imaginary part to Lab1

end;
end