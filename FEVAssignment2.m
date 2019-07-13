%% Loading Data and defining variables
[data] = xlsread('FEV Wesley.xlsx', 'Sheet1', 'A2:C343'); %taking till 18yrs

age = data(:,1);
fev = data(:,2);
smoke = data(:,3);
len = length(data);


%% Part (a)

j = 1;
k = 1;

for i = 1:len
    if smoke(i,:) == 1
        agesmoke(j,:) = age(i,:);
        fevsmoke(j,:) = fev(i,:);
        j = j+1;
    else
        agenonsmoke(k,:) = age(i,:);
        fevnonsmoke(k,:) = fev(i,:);
        k = k+1;
    end
end

tablerows = char('fev', 'age');

fprintf('         Sample Statistics for Smokers         \n')
fprintf('_________________________________________________\n')
fprintf('               Mean     Std Dev          \n')
fprintf('             _______   _________         \n')
fprintf('%-10s    %4.3f      %4.3f\n', tablerows(1,:), mean(fevsmoke), std(fevsmoke))
fprintf('%-10s   %4.3f      %4.3f\n', tablerows(2,:), mean(agesmoke), std(agesmoke))


tablerows = char('fev', 'age');
fprintf('\n')
fprintf('         Sample Statistics for Non-Smokers         \n')
fprintf('_________________________________________________\n')
fprintf('               Mean     Std Dev          \n')
fprintf('             _______   _________         \n')
fprintf('%-10s    %4.3f      %4.3f\n', tablerows(1,:), mean(fevnonsmoke), std(fevnonsmoke))
fprintf('%-10s   %4.3f      %4.3f\n', tablerows(2,:), mean(agenonsmoke), std(agenonsmoke))

%% Part (b)

whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [age, smoke];

statsa = regstats(fev, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
tstats = statsa.tstat.t;
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-3));

rowlabels = char('intercept','age', 'smoke');
fprintf('\n')
fprintf('        Descriptive Statistics for Part (b)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error    t stats          \n')
fprintf('                   _______   ________     _______ \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f  %10.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:), tstats(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);


%% Part (c)

agesqr = age.^2;
whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [age agesqr smoke];

statsa = regstats(fev, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
tstats = statsa.tstat.t;
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-4));

rowlabels = char('intercept', 'age', 'agesqr', 'smoke');
fprintf('\n')
fprintf('        Descriptive Statistics for Part (c)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error    t stats          \n')
fprintf('                   _______   ________     _______ \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f  %10.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:), tstats(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);


%% Part (d)

ageintsmoke = age.*smoke;
agesqr = age.^2;

whichstats = {'beta', 'covb', 'yhat', 'r', 'rsquare', 'mse', 'tstat', 'fstat'};
X = [age agesqr smoke ageintsmoke];

statsa = regstats(fev, X, 'linear', whichstats);

beta1 = statsa.beta;
SE1 = sqrt(diag(statsa.covb));
tstats = statsa.tstat.t;
res = statsa.r;
sse = norm(res,2)^2;
std = sqrt(sse/(len-5));

rowlabels = char('intercept', 'age', 'agesqr', 'smoke', 'ageIntsmoke');
fprintf('\n')
fprintf('        Descriptive Statistics for Part (d)         \n')
fprintf('_________________________________________________\n')
fprintf('                    Coeff    Std Error    t stats          \n')
fprintf('                   _______   ________     _______ \n')

for i = 1:length(beta1)
    fprintf('%-10s     %10.3f    %7.3f  %10.3f \n', rowlabels(i,:),...
        beta1(i,:), SE1(i,:), tstats(i,:));
end
fprintf('The R-square value is %0.3f \n', statsa.rsquare);
fprintf('The Standard Error is %0.3f \n', std);


%% Part (e)

% for Non-smoker, smoke = ageintsmoke = 0
Xnon = [ones(length(agenonsmoke),1) agenonsmoke agenonsmoke.^2 ...
    zeros(length(agenonsmoke),1) zeros(length(agenonsmoke),1) ];
plot(agenonsmoke, Xnon*beta1)
hold on
Xsmoke = [ones(length(agesmoke),1) agesmoke agesmoke.^2 ...
    ones(length(agesmoke),1) agesmoke ];
plot(agesmoke, Xsmoke*beta1, '--')
hold on 
scatter(age, fev,'g')
xlabel('Age');
ylabel('FEV');
legend('Non-Smoker','Smoker','Real Data');

%% Part (f)

% Hypothesis Testing

% Assuming aplha = 0.05
f_threshold = 3.022; % F(2,340)

% Restriction Parameters
X = [ones(len,1) age agesqr smoke ageintsmoke];
Res = [ 0 0 1 0 0 ; 0 0 0 0 1 ];
r = [ 0 ; 0 ];
k = 2;

f = (Res*beta1 - r)'*((Res*((X'*X)\Res'))\(Res*beta1 - r))/(k*statsa.mse);

if f < f_threshold
 fprintf('\n H_naught : beta_agesqr = beta_ageintsmoke =0 is accepted \n');
else
 fprintf('\n H_naught : beta_agesqr = beta_ageintsmoke = 0 is rejected \n');
end