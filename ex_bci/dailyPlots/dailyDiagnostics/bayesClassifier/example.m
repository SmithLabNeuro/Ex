clear; clc; close all;

load('fisherIrisData.mat');

[params,err] = crossvalidate_nb(X,y,'distType','sharedGauss','nFolds',10);

[pred_y,proj_y] = predict_nb(X,params);

fprintf('Cross-validated prediction accuracy: %.1f\n',(1-mean(err))*100)
fprintf('Training prediction accuracy: %.1f\n',100*sum(pred_y==y)./length(y))

