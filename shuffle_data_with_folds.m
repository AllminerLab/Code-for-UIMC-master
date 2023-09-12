function [datas,Yi,mapping] = shuffle_data_with_folds(raw_data,Y,folds,best_view )
% x is dv*n
n_views = length(raw_data); 
for i =1 :n_views
   raw_data{i}= raw_data{i}'; 
end
datas = {};
for a = 1 : n_views 
    n_view_a = size(raw_data{a}, 1);
    d_view_a = size(raw_data{a}, 2);
    mapping{a} = folds{a};  % random mapping

    datas{a} = zeros(size(raw_data{a}));
    
    for i = 1:n_view_a
       datas{a}(i,:) = raw_data{a}(mapping{a}(i),:);
    end
    

%     n_view_a = ceil(n_view_a - n_view_a * levels(a)); % remove samples after shuffle
%     datas{a} = datas{a}(1:n_view_a,:);
%     labels{a} = datas{a}(:,end);
%     datas{a} = datas{a}(1:n_view_a,1:end-1);
%     mappings{a} = mapping(1:n_view_a); % remove mappings
end
Yi = zeros(n_view_a,1);

    for i = 1:n_view_a
        Yi(i,:)  = Y(mapping{best_view}(i),:);
    end
for i =1 :n_views
   datas{i}= datas{i}'; 
end

end

