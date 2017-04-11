function fOut=p_ifft2d(fIn)

%**************************************************************************
%
%  fft2d computes a complex parallel 2D-FFT through 1-D Hartley Transform
%  using the row-column method.
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
%    Input, complex array fIn, the data to be transformed.
%
%    Output, complex array fOut, the transformed data.
%
%**************************************************************************
if (labindex==1) % Lab1 is master node
    [m,n]=size(fIn);
    br_row=floor(m/2); br_col=floor(n/2); % Points where array seperates
    
    labSend(fIn(br_row+1:end,:),2,0); % Send lowest rows to Lab2
    
    for i=1:br_row; % Trasform upper rows
        fIn(i,:)=cfht2fft(fIn(i,:),-1);
    end;
    
    labSend(fIn(1:br_row,br_col+1:end),2,1); % Send right upper block
    fIn(br_row+1:end,1:br_col)=labReceive(2,2); % Receive left down block
    
    for i=1:br_col; % Transform left columns
        fIn(:,i)=cfht2fft(fIn(:,i),1);
    end;
    
    tmp=labReceive(2,3); % Receive right columns from Lab2
    
    fOut=[fIn(:,1:br_col) tmp];
end;

if (labindex==2) % Lab2 executes from now on
    fIn=labReceive(1,0); % Receive rows for transform
    
    [m,n]=size(fIn);
    br_col=floor(n/2);
    
    for i=1:m; % Trasform lower rows
        fIn(i,:)=cfht2fft(fIn(i,:),-1);
    end;
    
    tmp=labReceive(1,1); % Receive right upper block
    labSend(fIn(:,1:br_col),1,2); % Send left down block
    
    fIn=[tmp;fIn(:,br_col+1:end)]; % Create columns for transform
    
    [~,n]=size(fIn); %Transform columns
    for i=1:n;
    	fIn(:,i)=cfht2fft(fIn(:,i),1);
    end;
    
    labSend(fIn,1,3); % Send right columns to Lab1
end;

end