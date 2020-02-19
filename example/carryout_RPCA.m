clear

addpath(genpath('./matlab_and_R_scripts')); 
tic
D=csvread('exampledata.csv',1,1);

%% there are no unobservable entries
%[A1,E1]= RPCA(D); % RPCA model

%% there exist unobservable entries
[m,n]=size(D);
omega=find(D~=3); % 3 represents missing
omegaC=find(D==3);
lambda=1/sqrt(max(m,n))*(1+3*length(omegaC)/(m*n));
[A1,E1]= extendedRPCA(D,omega,lambda); % the extended RPCA model

%% Integralization
AA1=int8(A1);
save('example_RobustClone.mat','AA1')
toc

