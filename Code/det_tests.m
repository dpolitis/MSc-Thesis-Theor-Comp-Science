
syms x y;

for maxPower=1:10
    fIn=create_bivariate_poly_matrix(maxPower); % Create test matrix

    fprintf('max power of polynomial matrix'); maxPower
    fprintf('\n');
    %t1=tic;
    %det(fIn); % Compute symbolic determinant
    %fprintf('symbolic determinant\n');
    %toc(t1)
    %fprintf('\n');
    
    t2=tic;
    p_interpolated_determinant(fIn,maxPower,maxPower);
    fprintf('interpolated determinant\n');
    toc(t2)
    fprintf('\n');
end;