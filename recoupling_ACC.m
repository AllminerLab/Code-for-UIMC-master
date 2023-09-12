function [ acc ] = recoupling_ACC(P,knn,gt_map)
% 1 构建epson近邻图以及记录index（best_view的)
% 2 P排序
% 3 P的排序后也就是tops是否属于对应的GT的KNN？
% 是的话count + 1 

[A,C] = sort(P,2,'descend');
[tops] = C(:,1);
cnt = 0;
for i = 1:size(P,1)   
    for j = 1:k
       if(map(1,i) == tops(i,j))
          cnt = cnt+1; 
       end
    end
end
acc = cnt/size(P,1);
end

