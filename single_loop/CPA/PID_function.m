function y = PID_function(k1,k2,k3,func_i,model_flag)
global m;
global nmi;
global gmi;

%% Benchmark function
if model_flag==1

    if func_i==1
        g_num=0.2;
        g_den=[1,-0.8];
        d=5; 
        T=1; 
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=1;
        n_den1=[1,-1];
        n_den2=[1,0.4];
        n_den=conv(n_den1,n_den2);
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);
    end
    
    if func_i==2
        g_num=0.08919;
        g_den=[1,-0.8669];
        d=12; 
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=0.08919;
        n_den=[1,-0.8669];
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);        
    end
    
    if func_i==3
        g_num=0.5108;
        g_den=[1,-0.9604];
        d=28;
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=0.5108;
        n_den=[1,-0.9604];
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d; 
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);
    end
    
    if func_i==4
        g_num=1;
        g_den=[1,-0.8];
        d=6;
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=[1,0.6];
        n_den1=conv([1,-0.5],[1,-0.6]);
        n_den2=[1,0.7];
        n_den=conv(n_den1,n_den2);
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);
    end
    
    if func_i==5
        g_num=1;
        g_den=[1,-0.8];
        d=6;
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=[1,-0.2];
        n_den1=conv([1,-1],[1,-0.3]);
        n_den2=conv([1,0.4],[1,-0.5]);
        n_den=conv(n_den1,n_den2);
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);   
    end
    
    if func_i==6
        g_num=1;
        g_den=[1,-0.8];
        d=6;
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=[1,0.6];
        n_den1=conv([1,-1],[1,-0.5]);
        n_den2=conv([1,-0.6],[1,0.7]);
        n_den=conv(n_den1,n_den2);
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);     
    end
    
    if func_i==7
        g_num=0.1;
        g_den=[1,-0.8];
        d=5;
        T=1; 
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=0.1;
        n_den1=conv([1,-1],[1,-0.3]);
        n_den2=[1,-0.6];
        n_den=conv(n_den1,n_den2);
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);       
    end
    
    if func_i==8
        g_num=0.1;
        g_den=[1,-0.8];
        d=3;
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=1;
        n_den=[1,-1];
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);
    end
    
    if func_i==9
        g_num=0.1;
        g_den=[1,-0.8];
        d=6;
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=0.1;
        n_den1=[1,-1];
        n_den2=[1,-0.7];
        n_den=conv(n_den1,n_den2);
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);       
    end
    
    if func_i==10
        g_num=0.1;
        g_den=[1,-0.8];
        d=3;
        T=1;
        g=tf(g_num,g_den,T,'Variable','z^-1','iodelay',d);
        n_num=sqrt(0.001);
        n_den1=[1,-1];
        n_den2=[1,0.2];
        n_den=conv(n_den1,n_den2);
        n=tf(n_num,n_den,T,'Variable','z^-1');
        m=8*d;
        nmi=impulse(n,m-1);
        gmi=step(g,m-1);       
    end
end

%% Output
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
Gmi=dh1\nmi;
y=Gmi'*Gmi;
end