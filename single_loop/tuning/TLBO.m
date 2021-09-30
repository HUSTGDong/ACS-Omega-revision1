function [Elite_antlion_fitness,var_y,var_u,teacher,G]=TLBO(ant_position,Np,Gm,Lowerbound,Upperbound,D,func_num,ruo)

func_init=PID_function(ant_position(1,1),ant_position(1,2),ant_position(1,3),func_num,ruo,1);

%***************Initialization**************%
G=1; 
ge=zeros(1,Gm);
Fit_Tr=zeros(Np,1);
St=ant_position;
St_new=zeros(Np,D);
Pop_Mean=zeros(Gm,D);
Pop_Variance=zeros(Gm,D);
Convergence_curve=zeros(1,200);

var_u_Np=zeros(Np,1);
var_y_Np=zeros(Np,1);

if size(Lowerbound,1) ==1 && size(Lowerbound,2)==1 %Check if the bounds are scalar
    Lowerbound=ones(1,D)*Lowerbound;
    Upperbound=ones(1,D)*Upperbound;
end
if size(Lowerbound,1) > size(Lowerbound,2) %Check if boundary vectors are horizontal or vertical
    Lowerbound=Lowerbound';
    Upperbound=Upperbound';
end

Elite_antlion_fitness=10;
ite_flag=9;
while abs(Elite_antlion_fitness-ite_flag)>10E-36
    %***************Teacher Phase**************%
    for i=1:D
        Sum_Column=0;
        for j=1:Np
            Sum_Column=Sum_Column+St(j,i);
        end
        Pop_Mean(G,i)=Sum_Column/Np;
    end
    %**********Calculate the variance of particles in each generation*********%
    for i=1:D
        Sum_Column_v=0;
        for j=1:Np
            Sum_Column_v=Sum_Column_v+(St(j,i)-Pop_Mean(G,i))^2;
        end
        Pop_Variance(G,i)=Sum_Column_v/Np;
    end
    for j=1:Np
        [Fit_Tr(j),var_y_Np(j),var_u_Np(j)]=PID_function(St(j,1),St(j,2),St(j,3),func_num,ruo,0);
    end
    [Y,Tridex]=min(Fit_Tr);
    Tr(G,:)=St(Tridex,:);
    for k=1:Np
        r1=rand;
        r2=rand;
        St_new(k,:)=St(k,:)+r1*(Tr(G,:)-((round(1+r2)*Pop_Mean(G,:))));
        [St_new]=Checkbound(St_new,Lowerbound,Upperbound,Np,D,G);
        if PID_function(St_new(k,1),St_new(k,2),St_new(k,3),func_num,ruo,0)<PID_function(St(k,1),St(k,2),St(k,3),func_num,ruo,0)
            St(k,:)=St_new(k,:);
        end
    end
    
    %***************Learner Phase**************%
    for i=1:Np
        dx=randperm(Np);
        if dx(1)==i 
            if PID_function(St(i,1),St(i,2),St(i,3),func_num,ruo,0)<PID_function(St(dx(2),1),St(dx(2),2),St(dx(2),3),func_num,ruo,0)
                St_new(i,:)=St(i,:)+rand*(St(i,:)-St(dx(2),:));
            else
                St_new(i,:)=St(i,:)+rand*(St(dx(2),:)-St(i,:));
            end
        else
            if PID_function(St(i,1),St(i,2),St(i,3),func_num,ruo,0)<PID_function(St(dx(1),1),St(dx(1),2),St(dx(1),3),func_num,ruo,0)
                St_new(i,:)=St(i,:)+rand*(St(i,:)-St(dx(1),:));
            else
                St_new(i,:)=St(i,:)+rand*(St(dx(1),:)-St(i,:));
            end
        end
        [St_new]=Checkbound(St_new,Lowerbound,Upperbound,Np,D,G);
        
        if PID_function(St_new(i,1),St_new(i,2),St_new(i,3),func_num,ruo,0)<PID_function(St(i,1),St(i,2),St(i,3),func_num,ruo,0)
            St(i,:)=St_new(i,:);
        end
    end
    
    LocalBest=Fit_Tr(Tridex);
    if G==1
        ge(1)=LocalBest;
    end
    if G>1
        if LocalBest<ge(G-1)
            ge(G)=LocalBest;
        else
            ge(G)=ge(G-1);
        end
    end
    
    var_y=var_y_Np(Tridex);
    var_u=var_u_Np(Tridex);
    Elite_antlion_fitness=LocalBest;
    % Update the convergence curve
    Convergence_curve(G)=Elite_antlion_fitness;
    
    if G>50
        ite_flag=Convergence_curve(G-50);
    end
          
    if mod(G,50)==0
        display(['At iteration ', num2str(G), ' the elite fitness is ', num2str(Elite_antlion_fitness)]);
    end
    
    G=G+1;
end
teacher=Tr(G-1,:);
end