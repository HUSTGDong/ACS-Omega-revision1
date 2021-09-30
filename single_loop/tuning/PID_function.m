function [J,var_y,var_u] = PID_function(k1,k2,k3,func_i,ruo,model_flag)
global m;
global nmi;
global gmi;
global nmi_step;

%% model
if model_flag==1
    if func_i==1
       g_num=[0.04133];
        g_den=[1,-0.8952];
        d=4; 
        T=1; 
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=[0.2];
        n_den=[1 -0.8952];
        n=tf(n_num,n_den,T,'Variable','z^-1');     
        m=12*d; 
        nmi=impulse(n,m*T-1);
        nmi_step=step(n,m*T-1);
        gmi=step(g,m*T-1);      
    end    
end
%% output variance
w=zeros(m,1);
w(1)=gmi(1)*k1;
w(2)=gmi(2)*k1+gmi(1)*k2;
for i=3:m
    w(i)=gmi(i)*k1+gmi(i-1)*k2+gmi(i-2)*k3;
end
W=zeros(m,m);
for j=1:m  
    for i=1:m
        if i>j
            W(i,j)=w(i-j+1);
        end
    end
end
I=eye(m);
dh1=I+W;

vmi=zeros(m,1);
vmi(1)=nmi_step(1)*k1;
vmi(2)=nmi_step(2)*k1+nmi_step(1)*k2;
for i=3:m
    vmi(i)=nmi_step(i)*k1+nmi_step(i-1)*k2+nmi_step(i-2)*k3;
end

Gmi=dh1\nmi_step;
Umi=dh1\vmi;
var_e=0.00001;
var_y=var_e*(Gmi'*Gmi);
var_u=var_e*(Umi'*Umi);
OV=var_y;
%% IAE
L=1000; 
y_set=3*ones(1,L);
e=zeros(1,L);
u=zeros(1,L);
y=zeros(1,L);

for i=1:L
    e(1)=y_set(1);
    u(1)=k1*e(1);
    if u(1)>16
        u(1)=16;
    end
    if u(1)<0
        u(1)=0;
    end
    
    e(2)=y_set(2)-y(1);
    u(2)=u(1)+k1*e(2)+k2*e(1);
    if u(2)>16
        u(2)=16;
    end
    if u(2)<0
        u(2)=0;
    end
    y(2)=0.8952*y(1);
    
    e(3)=y_set(3)-y(2);
    u(3)=u(2)+k1*e(3)+k2*e(2)+k3*e(1);
    if u(3)>16
        u(3)=16;
    end
    if u(3)<0
        u(3)=0;
    end
    y(3)=0.8952*y(2);
    
    e(4)=y_set(4)-y(3);
    u(4)=u(3)+k1*e(4)+k2*e(3)+k3*e(2);
    if u(4)>16
        u(4)=16;
    end
    if u(4)<0
        u(4)=0;
    end
    y(4)=0.8952*y(3);
    
    if i>4
        e(i)=y_set(i)-y(i-1);
        u(i)=u(i-1)+k1*e(i)+k2*e(i-1)+k3*e(i-2);
        
        if u(i)>16
            u(i)=16;
        end
        if u(i)<0
            u(i)=0;
        end
        
        y(i)=0.8952*y(i-1)+0.04133*u(i-4);
    end
end
IAE=sum(abs(e));
%% objective function
J=IAE+ruo*OV;
end