clear memory
clear;
addpath('./data'); 
addpath([pwd, '/funs']);
addpath('./average')
%% option
filename = 'yale_newdouble.mat';
% 'yale_newdouble.mat'  \ 'ORL.mat' \ 'MNIST.mat'

repNum = 1;  % Number of iterations
percentDel = 0.1;  % Missing percent

%% data_preprocessing
load(filename);
if strcmp('yale_newdouble.mat',filename)
    Y = gt; 
    best_view = 1;
elseif strcmp('ORL.mat',filename)
    Y=gt;
    best_view = 1;
elseif strcmp('MNIST.mat',filename)
    for i=1:length(X)
        X{i}=X{i}';
    end
    Y=truth;
    best_view = 3;
 end       

%% data pre processing
view_num = length(X);
t1=clock;
for v = 1:view_num
    X{v}=NormalizeData(X{v});
end

%% generate incomplete data
incomplete_path = ['./average/avg_incomplete/',filename(1:length(filename)-4),'_incomplete_per_', num2str(percentDel),'_best_',num2str(best_view) ,'.mat'];
load(incomplete_path)

%% shuffle data
uncoupled_path = ['./average/avg_uncoupled/',filename(1:length(filename)-4),'_uncoupled','_best_',num2str(best_view) , '.mat'];
load(uncoupled_path)

%% iteration
for iteration=1:repNum
    f = fix((iteration-1)/5 )+1;
    si_fold = shuffle_incomplete_fold(incomplete_folds{f},uncoupled_folds{f});
    [ISX,Y,mapping] = shuffle_data_with_folds(X,Y,uncoupled_folds{f},best_view); % shuffle
    [ISX, ~, W, Ne] = construct_incomplete(ISX,si_fold); % incomplete
  
  %% run
    [result,a1{iteration},a2{iteration},a3{iteration}] = UIMC(ISX,Y,W,filename,iteration,Ne,percentDel,best_view,mapping);
  %% result
    tempACC(iteration) = result(1);
    tempNMI(iteration) = result(2);
    tempPurity(iteration) = result(3);
    tempARI(iteration) = result(4);
    tempFscore(iteration) = result(5);
    tempPrecision(iteration) = result(6);
    tempRecall(iteration) = result(7);
end

ACC = [mean(tempACC),std(tempACC)];
NMI = [mean(tempNMI),std(tempNMI)];
Purity = [mean(tempPurity),std(tempPurity)];
ARI = [mean(tempARI),std(tempARI)];
Fscore = [mean(tempFscore),std(tempFscore)];
Precision = [mean(tempPrecision),std(tempPrecision)];
Recall = [mean(tempRecall),std(tempRecall)];

t2 = clock;
Time = etime(t2,t1)/repNum;
fprintf('\n----------------------------finished ------------------------\n'); 

disp(['filename = ', filename, ' best_view ',num2str(best_view),' percentDel ',num2str(percentDel)])
ACC = roundn(ACC,-4)
NMI = roundn(NMI,-4)
Purity = roundn(Purity,-4)
ARI = roundn(ARI,-4)
Fscore =roundn(Fscore,-4)
Precision =roundn(Precision,-4)
Recall =roundn(Recall,-4)
Time = roundn(Time,-4)

clear best_view enn knn f i incomplete_folds incomplete_path ISX iteration knn mapping Ne num_samples
clear num_Views percentDel repNum result si_fold t1 t2 uncoupled_folds uncoupled_path v view_num W X Y Z_ini
 clear tempPurity tempRecall Time tempACC tempARI tempFscore tempNMI tempPrecision 




