function [X,Z_ini,W,Ne] = construct_incomplete(X,ind_folds)
for iv = 1:length(X)
    X1 = X{iv};
%     X1 = NormalizeFea(X1,0);
%     X1 = NormalizeData(X{v});
    ind_0 = find(ind_folds(:,iv) == 0);
    ind_1 = find(ind_folds(:,iv) == 1);
    X1(:,ind_0) = 0;    % 缺失视角补0
    Y{iv} = X1;         % 一列一个样本
    % ------------- 构造缺失视角的索引矩阵 ----------- %
    linshi_W = eye(size(X1,2));
    linshi_W(:,ind_1) = [];
    W{iv} = linshi_W;
    Ne(iv) = length(ind_0); 
    % ---------- 初始KNN图构建 ----------- %
    X1(:,ind_0) = [];
    options = [];
    options.NeighborMode = 'KNN';
    options.k = 11;
    options.WeightMode = 'Binary';      % Binary  HeatKernel
    Z1 = full(constructW(X1',options));

    linshi_W = diag(ind_folds(:,iv));
    linshi_W(:,ind_0) = [];
    Z_ini{iv} = linshi_W*max(Z1,Z1')*linshi_W';
    clear Z1 linshi_W
end
clear X X1 ind_0
X = Y;
clear Y
% result = [X Z_ini W Ne];