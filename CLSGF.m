function [result] = MLAN_Clustering1(X,options)
% Multi-view Clustering with Consistent Local Structure-guided Graph Fusion
%% =====================  parameter =====================
c=options.c;    % the number of clusters
groundtruth=options.label;
k_tem = options.k_tem;  % the neighbor number of different views 
k = options.k;          % the neighbor number of the fused graph
for i=1:length(X)
    X{i}=full(X{i}');
    X{i}=double(X{i});
end
P= size(X,2); % view number
n = size(X{1},1);  % sample number
gamma = options.gamma;
NITER = 30;
%% =====================  Initialization =====================
W=ones(P,n)*(1/P);  % initialize W
A_p=cell(1,P); A_sum = zeros(n); 
for i = 1:P
    distX_initial(:,:,i) =  L2_distance_1( X{i}',X{i}' ) ; 
    [distXs_tem, idx_tem] = sort(distX_initial(:,:,i),2);
    A_p{i} = zeros(n);
    for ii = 1:n
        di = distXs_tem(i,2:k_tem+2);
        id = idx_tem(ii,2:k_tem+2);
        A_p{i}(ii,id) = (di(k_tem+1)-di)/(k_tem*di(k_tem+1)-sum(di(1:k_tem))+eps); % initialize similarity matrices A_p of different views
    end;
    A_sum = A_sum + diag(W(i,:).^2) *  A_p{i};
end
[dump, idx_a] = sort(-A_sum,2); 
S = zeros(n);
rr = zeros(n,1);
for i = 1:n
    di = -dump(i,2:k+2);
    rr(i) = sum(di(1:k)) - k*di(k);
    id = idx_a(i,2:k+2);
    S(i,id) = (di)/(sum(di(1:k))-k*di(k+1)+eps);               %initialize S
end;
lambda = mean(rr);   % setting lambda according to Eq. (28) in the paper
S = (S+S')/2;                                                
D = diag(sum(S));
L = D - S;
[V, ~, ~]=eig1(L, c, 0);    % initialize V
VV = L2_distance_1(V',V');
%% =====================  updating =====================
for iter = 1:NITER
    S = zeros(n);
    for i=1:n       
        idxa0 = idx_a(i,2:k+2);   
        dvi = VV(i,idxa0);   
        dxi = A_sum(i,idxa0);   
        gi =  2*dxi-0.5*gamma*dvi;
        Lambda_i = 1/k - 1/(2*k*sum(W(:,i).^2)+lambda)*sum(gi(1:k)); % setting ¦« according to Eq. (13) in the paper
        S(i,idxa0(1:k)) = gi(1:k)/(2*(sum(W(:,i).^2)+lambda)) + Lambda_i;   % updating S according to Eq. (11) in the paper     
    end
    S=max(S,0);
    
    h=[];
    for i = 1 : P
        h=[h ; sum( (A_p{i}-S).^2 ,2)'];
    end    
    for j=1:n
        tem=h(:,j).^(-1);
        W(:,j)=tem./sum(tem);  % update W according to Eq. (21) in the paper              
    end
    
    A_sum = zeros(n);  
    for i = 1 : P
        A_sum = A_sum + diag(W(i,:).^2) *A_p{i};
    end
     
    S = (S+S')/2;                                                      
    D = diag(sum(S));
    L = D-S;  
    V_old = V;
    [V, ~, ev]=eig1(L, c, 0);     % update V according to Eq. (22) in the paper    
    
    VV = L2_distance_1(V',V');
    [~, idx_a] = sort(-2*A_sum+0.5*gamma*VV,2);  
    % ------ gamma -----
     evs(:,iter+1) = ev;    
    thre = 1*10^-10;
    fn1 = sum(ev(1:c));      
    fn2 = sum(ev(1:c+1));
    if fn1 > thre
        gamma = 2*gamma;
    elseif fn2 < thre
        gamma = gamma/2;  V = V_old;
    else
        break;
    end;
end;

%% =====================  result =====================
[clusternum, y1]=graphconncomp(sparse(S)); y1 = y1';
if clusternum ~= c
    sprintf('Can not find the correct cluster number: %d', cluster_num)
end;
result = EvaluationMetrics(groundtruth', y1)*100;

fprintf('acc     nmi     pur     ar      f_sc\n');
fprintf('%0.2f   %0.2f   %0.2f   %0.2f   %0.2f\n', result(1), result(2), result(3),result(4),result(5));