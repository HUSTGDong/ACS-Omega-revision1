clear all 
clc
format long
%% Initialization
Function_Num=10;
dim=3;
SearchAgents_no=20; % Number of search agents
Max_iteration=50; % Maximum numbef of iterations
Method_Num=1;
Opera_Num=30;   

%% Process data
Time=zeros(Method_Num,Function_Num,Opera_Num);
Best_score=zeros(Method_Num,Function_Num,Opera_Num);
cg_curve=zeros(Method_Num,Function_Num,Opera_Num,200);
PID_para=zeros(Method_Num,Function_Num,Opera_Num,dim);
iter_num=zeros(Method_Num,Function_Num,Opera_Num);
ant_position=zeros(SearchAgents_no,dim);

%% Multigroup experiment
Method_num=1;
for i=1:Function_Num
    
    lb=-50;
    ub=50;
    
    for j=1:Opera_Num
       
        ant_position=initialization(SearchAgents_no,dim,ub,lb);

        tic;
        display(['problem',num2str(i),'  TLBO  ', num2str(j)]);
        [Best_score(Method_num,i,j),PID_para(Method_num,i,j,:),cg_curve(Method_num,i,j,:),iter_num(Method_num,i,j)]=TLBO(ant_position,SearchAgents_no,Max_iteration,lb,ub,dim,i);% 对蚂蚁和蚁狮初始值DOL
        time=toc
        Time(Method_num,i,j)=time;
        Method_num=Method_num+1;              
        
        if Method_num==Method_Num+1
            Method_num=1;
        end
        
    end
end

%% Data statistics
max_best=zeros(Method_Num,Function_Num);
min_best=zeros(Method_Num,Function_Num);
mean_best_score=zeros(Method_Num,Function_Num);
mean_time=zeros(Method_Num,Function_Num);
std_best_score=zeros(Method_Num,Function_Num);
for i=1:Method_Num
    for j=1:Function_Num
        mean_best_score(i,j)=mean(Best_score(i,j,:));
        mean_time(i,j)=mean(Time(i,j,:));
    std_best_score(i,j)=std(Best_score(i,j,:));
        max_best(i,j)=max(Best_score(i,j,:));
        min_best(i,j)=min(Best_score(i,j,:));
    end
end

%% Data saving
save 20210925_50_20

%% Mean and std of PID parameters
mean_PID = zeros(Method_Num,Function_Num,dim);
std_PID = zeros(Method_Num,Function_Num,dim);
for i=1:Method_Num
    for j=1:Function_Num   
        for k=1:dim
            mean_PID(i,j,k)=mean(PID_para(i,j,:,k));
            std_PID(i,j,k)=std(PID_para(i,j,:,k));
        end
    end  
end

mean_p=zeros(10,1);
mean_i=zeros(10,1);
mean_d=zeros(10,1);
std_p=zeros(10,1);
std_i=zeros(10,1);
std_d=zeros(10,1);
for i=1:Function_Num
    mean_p(i)=mean_PID(1,i,1);
    mean_i(i)=mean_PID(1,i,2);
    mean_d(i)=mean_PID(1,i,3);
    std_p(i)=std_PID(1,i,1);
    std_i(i)=std_PID(1,i,2);
    std_d(i)=std_PID(1,i,3);
end