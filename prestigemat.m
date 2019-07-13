[data] = xlsread('Prestige.xlsx', 'Sheet1');


%% Model 1
education  = data(:,1);
income = data(:,2);
women = data(:,3);
prestige = data(:,4);

educentered = education - mean(education);
wocentered = women - mean(women);
prescentered = prestige - mean(prestige);

%Linear Regression
%Model: y = b1 + x2b2 + x3b3 + x4b4

whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [educentered, prescentered, wocentered];

stats1 = regstats(income, X, 'linear', whichstats);

beta1 = stats1.beta;
SE1 = sqrt(diag(stats1.covb));
tstats1 = stats1.tstat.t;

rowlabels = char('intercept','educationcen', 'prestigecent', 'womencent');
fprintf('                        Descriptive Statistics         \n')
fprintf('_________________________________________________\n')
fprintf('                   coeff         SE       tstat          \n')
fprintf('                  ________     ______    _______         \n')

for i =1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f    %6.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:), tstats1(i,:));
end
stats1.mse
stats1.rsquare

scatter(stats1.yhat, stats1.r)



%% Second Regression Dropping education

education  = data(:,1);
income = data(:,2);
women = data(:,3);
prestige = data(:,4);

educentered = education - mean(education);
wocentered = women - mean(women);
prescentered = prestige - mean(prestige);

%Linear Regression
%Model: y = b1 + x2b2 + x3b3 + x4b4

whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [prescentered, wocentered];

stats2 = regstats(income, X, 'linear', whichstats);

beta2 = stats2.beta;
SE2 = sqrt(diag(stats2.covb));
tstats2 = stats2.tstat.t;

rowlabels = char('intercept','educationcent', 'prestigecent', 'womencent');
fprintf('              Descriptive Statistics         \n')
fprintf('_________________________________________________\n')
fprintf('                   coeff       SE        tstat          \n')
fprintf('                 _________   ________   _______         \n')

for i =1:length(beta1)
    fprintf('%-10s     %8.3f    %7.3f    %6.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:), tstats1(i,:));
end


