function aic = computeAIC(Xresid, numTrials, k)
    rss = sum((Xresid).^2);
    aic = numTrials * log(rss / numTrials) + 2 * k;
end