function pMatrixOut=create_bivariate_poly_matrix(maxPower)

%**************************************************************************
%
%  create_bivariate_poly_matrix computes a random 8x8 bivariate polynomial
%  matrix.
%
%  Discussion:
%
%    Firstly random coeffs are generated and transformed into random poly-
%    nomials by function create_bivariate_poly. The function is called 64
%    times to fill an 8x8 matrix with bivariate polynomials.
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
%    MatLab Help
%
%  Parameters:
%
%    Input, integer maxPower, the maximum polynomial degree.
%
%    Output, sym bivariate polynomial array fOut.
%
%**************************************************************************

pMatrixOut=sym(zeros(8,8)); syms x y; % Initialize variables

for i=1:8
    for j=1:8 % Fill polynomial matrix
        pMatrixOut(i,j)=create_bivariate_poly(maxPower);
    end;
end;

end

function pOut=create_bivariate_poly(maxPower)

syms x y; % Initialize variables

u=randi(4,1,1); % How many terms will the monomial have

% Monomial coeffs and powers
p=[randi(4,u,1) randi(maxPower,u,1) randi(maxPower,u,1) randi(4,u,1)];

pOut=0;

for i=1:u % Add monomials to create the polynomial
    pOut=pOut+p(i,1)*x^p(i,2)*y^p(i,3)+p(i,4);
end;

end

