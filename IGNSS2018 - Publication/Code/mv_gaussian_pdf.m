function pdf_value = mv_gaussian_pdf(x,u,C)
% Description
%
% Syntax: pdf_value = mv_gaussian_pdf(x,u,C)
% x: multi-variables vector 1*N
% u: mean vector 1*N
% C: Covariance Matrix N*N
% k: Dimension
x=x';
u=u';
k=length(x);
pdf_value = (1/sqrt ( (2*pi) ^ k * det(C) )) * exp (-1/2 * ( x - u)' * inv(C) * (x - u) );
    
end