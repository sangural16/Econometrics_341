[data, head] = xlsread('GradAdmiss.xlsx');

admit = data(:,1);
gre  = data(:,2);
gpa = data(:,3);
rank = data(:,4);

n = length(data);

rank1 = (rank == 1);
rank2 = (rank == 2);
rank3 = (rank == 3);
rank4 = (rank == 4);

X1 = [ gre gpa rank2 rank3 rank4 ];
[ betaProbit, devProbit, statsProbit ] = glmfit(X1, [admit ones(n,1)], ...
    'binomial', 'link', 'probit');
cov_prob = sqrt(diag(statsProbit.covb));

[ betaLogit, devLogit, statsLogit ] = glmfit(X1, [admit ones(n,1)], ...
    'binomial', 'link', 'logit');
cov_log = sqrt(diag(statsLogit.covb));

%% Minimization using fminsearch

%set options
options = optimset('MaxFunEvals', 2000, 'MaxIter', 2000);
%intial guess
b0 = [-2 0 0.5 -0.2 -0.61 -.53]';

tic
[b, fval, exitflag, output] = fminsearch('ProbitLogLike', ...
    b0, options, admit, [ones(n,1) X1]);
toc


%% Predicted Probablity of Admn

meanGRE = mean(gre);
meanGPA = mean(gpa);

Xpred = [ones(n,1) repmat(meanGRE, n, 1) repmat(meanGPA, n, 1) rank2 ...
    rank3 rank4];

meanPredProbRank1 = normcdf(Xpred*[betaProbit(1:3);0;0;0]);

[admit normcdf([ones(n,1) X1]*betaProbit)]

%% Covariate Effect

% GRE score increases by 25

Xgre = [ones(n,1) gre+25 gpa rank2 rank3 rank4];

PredXgre = normcdf(Xgre*betaProbit);
PredOldgre = normcdf([ones(n,1) X1]*betaProbit);

% change in Covariate Effect
 greCE = PredXgre - PredOldgre;

