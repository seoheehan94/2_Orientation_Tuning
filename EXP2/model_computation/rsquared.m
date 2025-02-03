
function r2 = rsquared(Xresid, Xorig)
r2 = 1 - sum((Xresid).^2)/sum((Xorig - mean(Xorig)).^2);
end