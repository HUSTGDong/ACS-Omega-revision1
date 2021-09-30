function [J,var_y,var_deltau] = PID_function(k1,k2,k3,func_i,ruo,model_flag)
global m;
global nmi1; % n1
global nmi2; % n2
global nmi1_diff; % n1
global nmi2_diff; % n2
% global nmi_step;
global gmi1; 
global gmi2; 
global gmi_step1; 
global gmi_step2; 

global g_num1;
global g_den1;
global g_num2;
global g_den2;

%% model
if model_flag==1

    if func_i==1
        g_num1=[0.04292];
        g_den1=[1,-0.9575];
        d1=7; 
        g_num2=[-0.5314];
        g_den2=[1,-0.6023];
        d2=3;    
        T=1;
        g1=tf(g_num1,g_den1,T,'Variable','z^-1','iodelay',d1);
        g2=tf(g_num2,g_den2,T,'Variable','z^-1','iodelay',d2);
        
        g_diff_num1=[1,-1];
        g_diff_den1=[1];
        g_diff=tf(g_diff_num1,g_diff_den1,T,'Variable','z^-1');
        
        n_num1=[1];
        n_den1=[1,-0.9575]; 
        n_num2=[1];
        n_den2=[1,-0.6023]; 

        n1=tf(n_num1,n_den1,T,'Variable','z^-1');
        n2=tf(n_num2,n_den2,T,'Variable','z^-1');
        
        m=16*d1; 
        nmi1=impulse(n1,m-1);
        nmi1_diff=impulse(n1*g_diff,m-1);
        nmi2=impulse(n2,m-1);
        nmi2_diff=impulse(n2*g_diff,m-1);
        gmi1=impulse(g1,m-1);
        gmi2=impulse(g2,m-1);
        gmi_step1=step(g1,m-1);
        gmi_step2=step(g2,m-1);
    end 
end

%% OV
S1=zeros(m,m);
S2=zeros(m,m);
delta_S1=zeros(m,m);
delta_S2=zeros(m,m);
F=zeros(m,m);
I=eye(m);
Tl=tril(ones(m));

for j=1:m  
    for i=1:m
        if i>j
            S1(i,j)=gmi_step1(i-j+1);
            S2(i,j)=gmi_step2(i-j+1);
            delta_S1(i,j)=gmi1(i-j+1);
            delta_S2(i,j)=gmi2(i-j+1);
        end
        if i==j && i<m
            F(i+1,j)=1;
        end
    end
end

var_deltau=0;


var_a1=0.00005;
var_a2=0.0005; 

w=(I+k3*delta_S2)\I;
W=k1*k3*delta_S1*w*S2+k2*k3*delta_S1*w*F*S2;
phi1=(I+W)\nmi1;
phi2=(I+W)\(delta_S1*w*nmi2);
var_y=phi1'*phi1*var_a1+phi2'*phi2*var_a2+2*phi1'*phi2*sqrt(var_a1)*sqrt(var_a2);

OV=var_y;
%% IAE
L=20000; 
y_set=3*ones(1,L);
e1=zeros(1,L);
e2=zeros(1,L);
u1=zeros(1,L);
u2=zeros(1,L);
y1=zeros(1,L);
y2=zeros(1,L);
u_lb=-16;
u_ub=0;

e1(1)=y_set(1);
u1(1)=k1*e1(1);
e2(1)=u1(1);
u2(1)=k3*e2(1);
if u2(1)>u_ub
    u2(1)=u_ub;
end
if u2(1)<u_lb
    u2(1)=u_lb;
end

e1(2)=y_set(2)-y1(1);
u1(2)=u1(1)+k1*e1(2)+k2*e1(1);
e2(2)=u1(2)-y2(1);
u2(2)=k3*e2(2);
if u2(2)>u_ub
    u2(2)=u_ub;
end
if u2(2)<u_lb
    u2(2)=u_lb;
end
y2(2)=-g_den2(2)*y2(1);
y1(2)=-g_den1(2)*y1(1);

e1(3)=y_set(3)-y1(2);
u1(3)=u1(2)+k1*e1(3)+k2*e1(2);
e2(3)=u1(3)-y2(2);
u2(3)=k3*e2(3);
if u2(3)>u_ub
    u2(3)=u_ub;
end
if u2(3)<u_lb
    u2(3)=u_lb;
end
y2(3)=-g_den2(2)*y2(2);
y1(3)=-g_den1(2)*y1(2);

e1(4)=y_set(4)-y1(3);
u1(4)=u1(3)+k1*e1(4)+k2*e1(3);
e2(4)=u1(4)-y2(3);
u2(4)=k3*e2(4);
if u2(4)>u_ub
    u2(4)=u_ub;
end
if u2(4)<u_lb
    u2(4)=u_lb;
end
y2(4)=-g_den2(2)*y2(3)+g_num2*u2(1);
y1(4)=-g_den1(2)*y1(3);

e1(5)=y_set(5)-y1(4);
u1(5)=u1(4)+k1*e1(5)+k2*e1(4);
e2(5)=u1(5)-y2(4);
u2(5)=k3*e2(5);
if u2(5)>u_ub
    u2(5)=u_ub;
end
if u2(5)<u_lb
    u2(5)=u_lb;
end
y2(5)=-g_den2(2)*y2(4)+g_num2*u2(2);
y1(5)=-g_den1(2)*y1(4);

e1(6)=y_set(6)-y1(5);
u1(6)=u1(5)+k1*e1(6)+k2*e1(5);
e2(6)=u1(6)-y2(5);
u2(6)=k3*e2(6);
if u2(6)>u_ub
    u2(6)=u_ub;
end
if u2(6)<u_lb
    u2(6)=u_lb;
end
y2(6)=-g_den2(2)*y2(5)+g_num2*u2(3);
y1(6)=-g_den1(2)*y1(5);

e1(7)=y_set(7)-y1(6);
u1(7)=u1(6)+k1*e1(7)+k2*e1(6);
e2(7)=u1(7)-y2(6);
u2(7)=k3*e2(7);
if u2(7)>u_ub
    u2(7)=u_ub;
end
if u2(7)<u_lb
    u2(7)=u_lb;
end
y2(7)=-g_den2(2)*y2(6)+g_num2*u2(4);
y1(7)=-g_den1(2)*y1(6);
    
for i=8:L
    e1(i)=y_set(i)-y1(i-1);
    u1(i)=u1(i-1)+k1*e1(i)+k2*e1(i-1);
    e2(i)=u1(i)-y2(i-1);
    u2(i)=k3*e2(i);
    if u2(i)>u_ub
        u2(i)=u_ub;
    end
    if u2(i)<u_lb
        u2(i)=u_lb;
    end
    y2(i)=-g_den2(2)*y2(i-1)+g_num2*u2(i-3);
    y1(i)=-g_den1(2)*y1(i-1)+g_num1*y2(i-7);
end

IAE=sum(abs(e1));

%% objective function
J=IAE+ruo*OV;
end