% The code for paper "Multi-view Clustering with Consistent Local Structure-guided Graph Fusion", IEEE Trans. Emerging Topics in Computational Intelligence, 2024.
clear all;
clc;
addpath('multiview_data'); 
addpath('tools'); 

%% setting
t=3; % 1-4 for different dataset, you only need to set this value for reproduction. 1£º3sources, 2£º100leaves, 3£ºMSRCv1, 4£ºCaltech101_20
 % demo
dataStr1={'3sources','100leaves','MSRCv1','Caltech101_20'};
k_tem=[150 100 100 80];
norm=[1 2 3 1];
options.dataStr1=dataStr1{t}; 
load(char(options.dataStr1));  options.c=numel(unique(label));
data=data_norm(data,norm(t)); 
options.label=label; 
options.c=numel(unique(label));  % the number of clusters
options.k_tem = k_tem(t);        % the neighbor number of different views 150 100 80 100  
options.k=15;                    % the neighbor number of the fused graph
options.gamma=1;                 % default

%% run
CLSGF(data,options);


