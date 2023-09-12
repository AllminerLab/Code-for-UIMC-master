function  [si_fold] = shuffle_incomplete_fold(in_fold, un_fold)
% in_fold = incomplete fold
% un_fold = uncoupled fold
num_views = size(in_fold,2);
n = size(in_fold,1);
for v = 1:num_views
    tmp_si_fold = zeros(size(in_fold,1),1);
    for ni = 1:n
        tmp_si_fold(ni) = in_fold(un_fold{v}(ni),v);
        % shuffle incomplete fold
    end
    si_fold(:,v)=tmp_si_fold;
    clear tmp_si_fold
end
end