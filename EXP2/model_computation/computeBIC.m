function bic = computeBIC(Xresid, numTrials, k)
    rss = sum((Xresid).^2);
    bic = numTrials * log(rss / numTrials) + k * log(numTrials);

end