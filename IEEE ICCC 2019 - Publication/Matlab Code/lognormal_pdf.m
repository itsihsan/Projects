function pdf_value = lognormal_pdf(x,u,C)
% Description

% Syntax: pdf_value = mv_gaussian_pdf(x,u,C)
% x: multi-variables vector 1*N
% u: mean vector 1*N
% C: Variance Matrix
% k: Dimension

%x=x'; % 6 X 20000 ----- Y vector
%u=u'; % 6 X 20000 ----- U/V vector

N=length(x);

for i=1:1:N
    pdf_value(i,:) = exp((-(x(i,:)-u(i,:)).^2)/(2*C));
end
end


% % %  function pdf_value = lognormal_pdf(x,u,C)
% % % % Description
% % % %
% % % % Syntax: pdf_value = mv_gaussian_pdf(x,u,C)
% % % % x: multi-variables vector 1*N
% % % % u: mean vector 1*N
% % % % C: Covariance Matrix
% % %
% % % x=x'; % Y
% % % u=u'; % U or V
% % %
% % % % pdf of a LOG normal distribution -- some of the tersm deleted due to LRT
% % % % pdf_value = (1./det_1) * exp((-num_2.^2)/(2*C));
% % % pdf_value = exp((-(x-u).^2)/(2*C));
% % %
% % % end