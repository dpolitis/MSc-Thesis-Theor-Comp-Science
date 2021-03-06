function coeffs=p_interpolated_determinant(A,l,k)

%**************************************************************************
%
%  This is a test function to compute the determinant of 2-D polynomial
%  matrix using parallell interpolation and Fast Hartley Transform.
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
%**************************************************************************

if (labindex==1) % Lab1 is master node
    syms x y % Initialize variables
    
    [m n]=size(A); % Polynomial matrix dimensions
    
    M1=l*m; M2=k*n;
    
    R=(M1+1)*(M2+1); % Number of interpolation points
    
    W1=exp(2*pi*1i/(M1+1));
    W2=exp(2*pi*1i/(M2+1)); % General case where M1~=M2
    
    r1=0:M1; r2=0:M2; % Calculate Interpolation points
    u1=W1.^(-r1); u2=W2.^(-r2);
    
    detA=zeros(M1+1,M2+1); % Genarate indices
    i=1:M1+1; tmp=ones(1,M2+1); i=kron(i,tmp);
    j=1:M2+1; tmp=ones(1,M1+1); j=kron(tmp,j);
    tmp=[i;j];
    
    labSend(tmp,2,1);
    labSend(u1,2,2); labSend(u2,2,3);
    labSend(R,2,4);
    
    for k=1:floor(R/2) % Evaluate determinant on Interpolation points
        i=tmp(1,k); j=tmp(2,k);
        x=u1(i); y=u2(j);
        detA(i,j)=det(eval(A));
    end;
    
    temp=labReceive(2,5);
    detA=[detA;temp];
    
    coeffs=ifft2d(detA); % Compute polynomial coeficients
end;

if (labindex==2) % Lab2 goes from now on
    syms x y % Initialize variables
    
    tmp=labReceive(1,1);
    u1=labReceive(1,2); u2=labReceive(1,3);
    R=labReceive(1,4);
    
    for k=floor(R/2)+1:R % Evaluate determinant on Interpolation points
        i=tmp(1,k); j=tmp(2,k);
        x=u1(i); y=u2(j);
        detA(i,j)=det(eval(A));
    end;
    
    labSend(detA(floor(R/2)+1:end,:),1,5);
end;

end