function[My_result,a1,a2,a3,K_Z,P_M,k_acc]= UIMC(X,Y,W,filename,iteration,Ne,percentDel,best_view,mapping)
%% Initialization
gt = Y;
clear Y 
cls_num = length(unique(gt));
Label = double(gt);
num_views = length(X); 
N = size(X{1},2); 
P_M = [];
K_Z= [];
beta = 1 ./ num_views;
gamma = 1;
tmp_Z = rand(N,N);
tmp_P = rand(N,N);
tmp_M = rand(N,N);

for v = 1:num_views
    Z{v} = bsxfun(@times, tmp_Z, 1./sum(tmp_Z));
    Q{v} = zeros(N,N);
    E{v} = zeros(size(X{v},1),N);
    S{v} = zeros(size(X{v},1),Ne(v));
    W{v} = W{v}'; % W is the indicator
    P{v} = bsxfun(@times, tmp_P, 1./sum(tmp_P));
    M{v} = bsxfun(@times, tmp_M, 1./sum(tmp_M)); 
    K{v} = zeros(N,N);
    J{v} = zeros(size(X{v},1),N);
    G{v} = zeros(N,N);
    Y{v}=X{v}+S{v}*W{v};
    alpha{v} = 0;
end
I = eye(N,N);
sX = [N, N, num_views];
iter = 0;Isconverg = 0; num = 0;
epson = 1e-5;
mu = 10e-5; max_mu = 10e10; pho_mu = 2; 

%% iteration 
while(Isconverg ==0)
    iter = iter +1;
    fprintf('---------processing epoch %d ----- iter %d--------\n', iteration,iter); 
    num = num + 1;
    for v = 1:num_views
        if v == best_view
            alpha{v} = 0;
        else
            alpha{v} = 1 ./ (2 * sqrt(  norm( Z{v}*P{v}-Z{best_view},"fro")  ));
        end
    end
    for v = 1: num_views
        
    %% Update Z^k
        Y{v}=X{v}+S{v}*W{v};
        tmp_A = mu*(Y{v})'*(Y{v});
        tmp_B = 2*alpha{v}*M{v}*M{v}'+mu*P{v}*P{v}';
        tmp_C = 2*alpha{v}*Z{best_view}*M{v}' + mu*(Y{v})'*(Y{v}) - mu*(Y{v})'*E{v} + (Y{v})'*J{v} + mu*K{v}*P{v}' - Q{v}*P{v}';
        Z{v} =  lyap(tmp_A,tmp_B,-tmp_C);
        clear tmp_A tmp_B tmp_C tmp_Z
        
    %% Update M^v
        if v == best_view
            M{v}=eye(N,N);
        else
          M{v} = pinv(2*alpha{v}*Z{v}'*Z{v}+mu*I)*(2*alpha{v}*Z{v}'*Z{best_view}+mu*P{v}+G{v});
        end
        
    %% Update P^v
        if v == best_view
            P{v} = eye(N,N); 
        else
            P{v} = pinv(Z{v}'*Z{v}+I)*(M{v}+Z{v}'*K{v} - Z{v}'*Q{v}./mu - G{v}./mu);
        end
        
   %% S^v sub problem 
        tmp_WZ = W{v}-W{v}*Z{v};
        tmp_X = X{v}-X{v}*Z{v}-E{v}+J{v}./mu;
        tmp_S = - (mu*tmp_X*tmp_WZ')*pinv((2*beta*eye(Ne(v)) + mu * (tmp_WZ * tmp_WZ')));
        S1 = zeros(size(tmp_S));
        for is = 1:size(tmp_S,2)
           S1(:,is) = EProjSimplex_new(tmp_S(:,is));
        end
            S{v} = S1;
        clear tmp_S S1 tmp_WZ tmp_X
         
    %% Update Jv Gv
        tmpY{v} = (X{v}+S{v}*W{v});
        J{v} = J{v} + mu * (tmpY{v}-tmpY{v} * Z{v}-E{v});
        G{v} = G{v} + mu * (P{v}-M{v});
        tmp_PZ{v} = Z{v}*P{v};  
    end
    
    %% Update E^k
    F = [];
    for v = 1 : num_views
        F = [F; tmpY{v}-tmpY{v}*Z{v}+J{v}/mu];
    end
    [Econcat] = solve_l1l2(F,gamma/mu);
    for v = 1:num_views
        if v == 1
            concatlength = 0;
        end
         E{v} = Econcat((concatlength+1):(concatlength+size(tmpY{v},1)),:);  
         concatlength = concatlength + size(tmpY{v});
    end    
    
    %% Update K
    Z_tensor = cat(3, tmp_PZ{:,:});
    Q_tensor = cat(3,Q{:,:});
    z = Z_tensor(:);
    q = Q_tensor(:);
    [k, ~] = wshrinkObj_weight(z + 1/mu*q,1/mu,sX,0,3);
    K_tensor = reshape(k,sX);
    %% Update Q
    q = q+mu*(z-k); 
    
    %% Update mu
    mu = min(mu*pho_mu, max_mu);
    
    %% Converge Condition
    Q_tensor = reshape(q,sX);
    Isconverg = 1;
    max_K_Z = 0;
    
    for v = 1: num_views
        K{v}=K_tensor(:,:,v);
        Q{v}=Q_tensor(:,:,v);    
        
        if (norm(K{v}-Z{v}*P{v},inf)>epson)
            history.norm_K= norm(K{v}-Z{v}*P{v},inf);
            fprintf('norm_KZ_v%1.0f %7.10f\n', v,history.norm_K);
            Isconverg = 0;
            max_K_Z=max(max_K_Z,history.norm_K);
        end
        a1{v}(iter)= mappingsACC(M{v},mapping{best_view},1); 
        a2{v}(iter) = mappingsACC(M{v},mapping{best_view},3);
        a3{v}(iter) = mappingsACC(M{v},mapping{best_view},10);
    end
 
end
    tmp=zeros(N,N);
    for v=1: num_views
        tmp = tmp + tmp_PZ{v};
    end
    Clus = SpectralClustering(1/(2*num_views)*(abs(tmp)+abs(tmp')), cls_num );

    [result(num,:)] = ClusteringMeasure1(Label, Clus);
    My_result = result(num,:);   
    


   
    
    