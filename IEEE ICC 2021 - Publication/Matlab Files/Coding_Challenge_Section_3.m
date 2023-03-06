%%% Please note that claculation for N=20 will follow the same steps as for
%%% N=10. The time is close and my PC is taking really long time to
%%% finalize "factorial(20) = 2.4329e+18" calculations for N=20

N = 10; % Total coins
drawnValues = N; % Sample space
for i = 1:factorial(N)
    y = (randsample(N,drawnValues,false))';
    absDiff = abs(y(1:drawnValues-1) - y(2:drawnValues));
    totPayment(i) = y(1)+sum(absDiff);
    y=[];
    absDiff=[];
end
digitsOld = digits(10)
meanTotPayment = vpa(mean(totPayment))
stdTotPayment = vpa(std(totPayment))
nGreatThan45 = vpa(length(totPayment(totPayment>=45)));
probGreatThan45 = vpa(nGreatThan45/length(totPayment))