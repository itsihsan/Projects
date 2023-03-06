function pdf_value = gaussian_pdf(x,u,C)
% Description
% gaussian_pdf(Y_G,V,R);
% Syntax: pdf_value = mv_gaussian_pdf(x,u,C)
% x: multi-variables vector 1*N
% u: mean vector 1*N
% C: Covariance Matrix N*N
% k: Dimension

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % Lognroaml case C is a row matrix 1 X N.
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % Making variance covariance matrix.
% % % 
% % % before_bracket = [];
% % % var_cov = [];
% % % in_bracket = [];
% % % in_bracket = exp(diag(C))- eye(length(C))
% % % 
% % % rows     =  size(x,1) % Number of rows/measurments.
% % % columns  =  size(x,2) % Number of columns/verifiers.
% % % 
% % % for ii=1:rows
% % %     for jj=1:columns
% % %         Expect_value(ii,jj) = exp(u(ii,jj) + (0.5*C(1,jj)))
% % %     end
% % %     before_bracket = Expect_value(ii,:) * (Expect_value(ii,:))'
% % %     var_cov = before_bracket .* in_bracket
% % %     pdf_value(ii,:)= exp (-1/2 * ( x(ii,:) - u(ii,:)) * inv(var_cov) * (x(ii,:) - u(ii,:))' )
% % %     
% % %     before_bracket = [];
% % %     var_cov = [];
% % % end


%%%% Ziqing previous code
% C = diag(C) % command added for log normal case.


x=x'; % 6 X 20000
u=u'; % 6 X 20000
k=length(x);
N=length(x);
for i=1:1:N
    pdf_value(i,:)= exp (-1/2 * ( x(:,i) - u(:,i))' * inv(C) * (x(:,i) - u(:,i)) );
end
end