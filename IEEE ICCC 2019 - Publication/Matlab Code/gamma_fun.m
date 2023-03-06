function gam = gamma_fun(yy,dd)

RSS_ref = -65.7;
d_ref = 37;

for i=1:length(yy)
    gam(i) = (yy(i)-RSS_ref)/(10*log10(d_ref/dd(i)));
end