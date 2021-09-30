clear all 
clc
format long
%% initialization
Function_Num=1;
dim=3;
SearchAgents_no=40; % Number of search agents
Max_iteration=500; % Maximum numbef of iterations

%% tuning
ruo=[0,1000000,10000000,100000000];
ruo_num=length(ruo);
Time=zeros(Function_Num,ruo_num);
var_u=zeros(Function_Num,ruo_num);
var_y=zeros(Function_Num,ruo_num);
Best_score=zeros(Function_Num,ruo_num);
PID_para=zeros(Function_Num,ruo_num,dim);
iter_num=zeros(Function_Num,ruo_num);
ant_position=zeros(SearchAgents_no,dim);

for i=1:Function_Num
    
    lb = [-50,-50,-50];
    ub = [50,50,50];
    
    for j=1:ruo_num
        
        ant_position=initialization(SearchAgents_no,dim,ub,lb);
        antlion_position=initialization(SearchAgents_no,dim,ub,lb);

        tic;
        display(['problem',num2str(i),'  TLBO  ', num2str(j)]);
        [Best_score(i,j),var_y(i,j),var_u(i,j),PID_para(i,j,:),iter_num(i,j)]=TLBO(ant_position,SearchAgents_no,Max_iteration,lb,ub,dim,i,ruo(j));% 对蚂蚁和蚁狮初始值DOL
        time=toc
        Time(i,j)=time;      
    end
end

%% data saving
save 20210927