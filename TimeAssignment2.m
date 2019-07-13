%% Loading Data and defining variables
[data] = xlsread('time.xlsx', 'Sheet1');

y1 = data(:,2);
y2 = data(:,3);
y3 = data(:,4);
x2 = data(:,5);
x3 = data(:,6);

len = length(data);
Dt = zeros(len,1);

for i = 1:20
    if data(i,1) >= 1939 && data(i,1) <= 1945
        Dt(i) = 1;
    else
        Dt(i) = 0;
    end
end

Ct = 1-Dt;


%% Part (a)
whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [x2, x3];

statsa = regstats(y1, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-3));

rowlabels = char('intercept','x2', 'x3');
fprintf('        Descriptive Statistics for Part (a)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error          \n')
fprintf('                   _______   ________            \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);


%% Part (b)
whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [x2, x3];

statsa = regstats(y3, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-3));

rowlabels = char('intercept','x2', 'x3');
fprintf('\n');
fprintf('        Descriptive Statistics for Part (b)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error          \n')
fprintf('                   _______   ________            \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);


%% Part (c)
whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [Ct, Dt, x2, x3];

statsa=regstats(y1,X,[1 0 0 0 ; 0 1 0 0 ;0 0 1 0 ;0 0 0 1], whichstats);
cov_var = statsa.covb; %% variance-covariance matrix
beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-4));

cov_var_estimate=(statsa.mse)*(pinv(X'*X));
fprintf('\n');
fprintf('The estimated Var-Covar matrix:\n');
disp(cov_var_estimate);

rowlabels = char('Ct', 'Dt', 'x2', 'x3');
fprintf('\n');
fprintf('        Descriptive Statistics for Part (c)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error          \n')
fprintf('                   _______   ________            \n')
for i = 1: length(beta1)
fprintf('%-10s     %10.3f    %7.3f \n', rowlabels(i,:),...
    beta1(i), SE1(i));
end

fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);

%% Part (d)
covar = cov_var_estimate;
var_delta = covar(1,1) + covar(2,2) - 2*covar(1,2);
fprintf('\n');
fprintf('Estimated var(1,1)+Estimated var(2,2)-2*Estimated var(2,1)=%7.3f\n',...
    var_delta);

%% Part (e)

whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [Dt x2 x3];

statsa = regstats(y1, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-4));

rowlabels = char('intercept', 'Dt', 'x2', 'x3');
fprintf('\n');
fprintf('        Descriptive Statistics for Part (e)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error          \n')
fprintf('                   _______   ________            \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);


%% Part (f)
x4 = x3.*Dt;
x5 = x3.*Ct;

whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [Dt x2 x3 x4];

statsa = regstats(y2, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-5));

rowlabels = char('intercept', 'Dt', 'x2', 'x3', 'x4');
fprintf('\n');
fprintf('        Descriptive Statistics for Part (f)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error          \n')
fprintf('                   _______   ________            \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);


%% Part (g)
whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [Ct Dt x2 x4 x5];

statsa = regstats(y2, X, [1 0 0 0 0; 0 1 0 0 0; 0 0 1 0 0; 0 0 0 1 0; 0 0 0 0 1], whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-5));

cov_var_estimate=(statsa.mse)*(pinv(X'*X));
fprintf('\n');
fprintf('The estimated Var-Covar matrix:\n');
disp(cov_var_estimate);

rowlabels = char('Ct', 'Dt', 'x2', 'x4', 'x5');
fprintf('\n');
fprintf('        Descriptive Statistics for Part (g)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error          \n')
fprintf('                   _______   ________            \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);

% Hypothesis Testing

% Assuming aplha = 0.05
% F(2,15)
f_threshold = 3.6823;
% Restriction Parameters
Res = [ 0 0 0 1 0 ; 0 0 0 0 1 ];
r = [ 0 ; 0 ];
k = 2;

f = (Res*beta1 - r)'*((Res*((X'*X)\Res'))\(Res*beta1 - r))/(k*statsa.mse);

if f < f_threshold
 fprintf('\n H_naught : beta4 = beta5 = 0 cannot be rejected \n');
else
 fprintf('\n H_naught : beta4 = beta5 = 0 is rejected \n');
end



%% Part (h)
x6 = x2.*Dt;

whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [Dt x2 x3 x4 x6];

statsa = regstats(y3, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-6));

cov_var_estimate=(statsa.mse)*(pinv(X'*X));
fprintf('\n');
fprintf('The estimated Var-Covar matrix:\n');
disp(cov_var_estimate);

rowlabels = char('intercept', 'Dt', 'x2','x3', 'x4', 'x6');
fprintf('\n');
fprintf('        Descriptive Statistics for Part (h)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error          \n')
fprintf('                   _______   ________            \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);