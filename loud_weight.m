function R = loud_weight(f)

R = 12200^2.*f.^4/(f.^2+20.6^2)./sqrt((f.^2+107.7^2).*(f.^2+737.9^2))./(f.^2+12200^2);
