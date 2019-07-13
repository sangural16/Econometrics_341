%%Simulation Example
%Model +: y = b1 + b2x2 + b3x3 + eps

n = 1000;
b1 = 5; b2 = 3; b3 = 2;
x2 = normrnd(0,1,n,1);
x3 = normrnd(0,1,n,1);
eps = normrnd(0,1,n,1);

y = b1.*ones(n,1) + b2.*x2 + b3.*x3 + eps;
[y x2 x3];

%estimating the model
whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [x2 x3];
stats2 = regstats(y,X,'linear', whichstats);

disp('    beta est    std err    t-stat');
disp([stats2.beta sqrt(diag(stats2.covb)) stats2.tstatans
    .t])



