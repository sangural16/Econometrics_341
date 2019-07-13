%% Whites Heteroskedastic Test

covOmega = stats2.r*stats2.r';
Omega = diag(diag(covOmega));

invOmega = Omega\eye(n);

X1 = [ones(n,1) X];

bfgls = (X1'*invOmega*X1)\(X1'*invOmega*y);

minx2 = min(x2);
maxx2 = max(x2);

x2seq = minx2:.05:maxx2;
m = length(x2seq);
Xpred = [ones(m,1) x2seq' mean(x3).*ones(m,1)];
betahat = stats2.beta;
ypred = Xpred*betahat;

scatter(x2,y)
hold on
plot(x2seq, ypred, 'k')