function pdf_value = gaussian_pdf(x,u,C)
% Description
% gaussian_pdf(Y_G,V,R);
% Syntax: pdf_value = mv_gaussian_pdf(x,u,C)
% x: multi-variables vector 1*N
% u: mean vector 1*N
% C: Covariance Matrix N*N
% k: Dimension

% C = diag(C) % command added for log normal case.


x=x'; % 6 X 20000
u=u'; % 6 X 20000
k=length(x);
N=length(x);
for i=1:1:N
    pdf_value(i,:)= exp (-1/2 * ( x(:,i) - u(:,i))' * inv(C) * (x(:,i) - u(:,i)) );
end
end