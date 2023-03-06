function  [RSS]= add_shadowing (u,R)

% Calculate the correlated shadowing and add it to the RSS mean value.
RSS_mean = u;
num = length(u);
shadowing_variance = 5;  % Shadowing Variance [Sigma]
correlated_shadowing = shadowing_variance * mvnrnd ( zeros (1,num) , R );
RSS = RSS_mean + correlated_shadowing;

