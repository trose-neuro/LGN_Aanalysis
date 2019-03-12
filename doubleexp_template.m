function [a b c d gof]=doubleexp_dec(time, intensity,fo_, st_)

xfit=time;
yfit=intensity;

set(fo_,'Startpoint',st_);  set(fo_,'MaxFunEvals',600); set(fo_,'MaxIter',400)


ft_ = fittype('a*exp(-b*x) + c*exp(-d*x)',...
    'dependent',{'y'}, 'independent',{'x'},...
    'coefficients', {'a', 'b', 'c', 'd'},  'options', fo_);

% Fit this model using new data
[cf_ gof]= fit(xfit,yfit,ft_ ,'Startpoint',st_);
rsquare_bleach = gof.rsquare;

coeffs = coeffvalues(cf_);
a = (coeffs(1));
b = (coeffs(2));
c = (coeffs(3));
d = (coeffs(4));
end
