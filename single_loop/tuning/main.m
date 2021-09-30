clear all 
clc
format long
%% Initialization
Function_Num=1;
dim=3;
SearchAgents_no=20; % Number of search agents
Max_iteration=500; % Maximum numbef of iterations

%% PID parameters
ruo=[0,100000,250000,1000000];
ruo_num=length(ruo);
Time=zeros(Function_Num,ruo_num);
var_u=zeros(Function_Num,ruo_num);
var_y=zeros(Function_Num,ruo_num);
Best_score=zeros(Function_Num,ruo_num);
PID_para=zeros(Function_Num,ruo_num,dim);
iter_num=zeros(Function_Num,ruo_num);
ant_position=zeros(SearchAgents_no,dim);

i=1;

lb = [-50,-50,-50];
ub = [50,50,50];

for j=1:ruo_num
    
    ant_position=initialization(SearchAgents_no,dim,ub,lb);
    
    tic;
    display(['problem',num2str(i),'  TLBO  ', num2str(j)]);
    [Best_score(i,j),var_y(i,j),var_u(i,j),PID_para(i,j,:),iter_num(i,j)]=TLBO(ant_position,SearchAgents_no,Max_iteration,lb,ub,dim,i,ruo(j));
    time=toc
    Time(i,j)=time;
    
end

%% Data saving
save 20210926