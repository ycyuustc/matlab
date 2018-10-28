function [result]=pl_gaudin_insert_epsilon(h,mu,Temperature)

result=zeros(1,9);
global T;
T=Temperature;

global v_k1;
global v_k2;
global v_ep1;
global v_ep2;
global v_ep1_dmu;
global v_ep2_dmu;
global v_ep1_dh;
global v_ep2_dh;
global v_ep1_dmudmu;
global v_ep2_dmudmu;
global v_ep1_dhdh;
global v_ep2_dhdh;
global v_ep1_dT;
global v_ep2_dT;
global v_ep1_dTdT;
global v_ep2_dTdT;
global v_ep1_dTdmu;
global v_ep2_dTdmu;
global v_ep1_dTdh;
global v_ep2_dTdh;

 global v_f_ep1;
 global v_f_ep2;
global v_g_ep1;
global v_g_ep2;
global v_h_ep1;
global v_h_ep2;
global v_h1_ep1;
global v_h1_ep2;
global v_h2_ep1;
global v_h2_ep2;


a1=@(x) 1./(2*pi)./(1/4+x.^2);
a2=@(x) 2./(2*pi)./(1+x.^2);
Num_k=1000;

Num_k1=Num_k;
kc1=0.5;

v_k1=linspace(0,kc1,Num_k1);
v_length=v_k1(2:Num_k1)-v_k1(1:Num_k1-1);
v_k1_nt=zeros(1,Num_k1);
v_k1_nt(1:Num_k1-1)=1/2*v_length;
v_k1_nt(2:Num_k1)=v_k1_nt(2:Num_k1)+1/2*v_length;

Num_k2=Num_k-100;
kc2=0.5;

v_k2=linspace(0,kc2,Num_k2);
v_length=v_k2(2:Num_k2)-v_k2(1:Num_k2-1);
v_k2_nt=zeros(1,Num_k2);
v_k2_nt(1:Num_k2-1)=1/2*v_length;
v_k2_nt(2:Num_k2)=v_k2_nt(2:Num_k2)+1/2*v_length;


[tmk,tmq]=meshgrid(v_k2,v_k2);
m_a2_22=a2(tmk-tmq)+a2(tmk+tmq);

[tmk,tmq]=meshgrid(v_k1,v_k2);
m_a1_21=a1(tmk-tmq)+a1(tmk+tmq);  %  2 transform to 1

[tmk,tmq]=meshgrid(v_k2,v_k1);
m_a1_12=a1(tmk-tmq)+a1(tmk+tmq); %  1 transfrom to 2

v_ep1_ori=v_k1.^2-mu-h/2;
v_ep2_ori=2*v_k2.^2-2*mu-1/2;

v_ep1=v_ep1_ori;
v_ep2=v_ep2_ori-v_k1_nt.*fn_F(v_ep1)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_pre=v_ep1;
    v_ep2_pre=v_ep2;
    
    v_ep1=v_ep1_ori-v_k2_nt.*fn_F(v_ep2)*m_a1_21;
    v_ep2=v_ep2_ori-v_k1_nt.*fn_F(v_ep1)*m_a1_12 ...
        -v_k2_nt.*fn_F(v_ep2)*m_a2_22;
    
    flag1=max(abs(v_ep1-v_ep1_pre));
    flag2=max(abs(v_ep2-v_ep2_pre));
    flag=max(flag1,flag2);
    
end

result(1)=-1/pi/2*sum(v_k1_nt.*fn_F(v_ep1))...
    -1/pi*sum(v_k2_nt.*fn_F(v_ep2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   insert point 
NUM_int=41;

bound=log(1000/(T^2));                         % ·ÉÑ©ÖÐÕªÒ»Æ¬ÂäÒ¶
tv_add=linspace(-bound*0.6,bound*0.6,NUM_int);

%v_int_point=T*[-10,-4.80,-3,-1.32,-0.35,0,0.35,1.32,3,4,80,10];
v_int_point=T*tv_add;

v_k1_ins=zeros(1,NUM_int);
v_ep1_ins=zeros(1,NUM_int);
v_k2_ins=zeros(1,NUM_int);
v_ep2_ins=zeros(1,NUM_int);
length_vk1_ins=0;
length_vk2_ins=0;

for i_ins=1:NUM_int
    
    find_e=v_int_point(i_ins);  
    tvq=v_k1;
    tve=v_ep1;
    
    v_index=find((tve(1:end-1)-find_e).*(tve(2:end)-find_e)<0);
    
    for i2=1:length(v_index)
        
        t_index=v_index(i2)+1;
        
        q1=tvq(t_index-1);
        q2=tvq(t_index);
        
        e1=tve(t_index-1);
        e2=tve(t_index);
        
        e1=e1-find_e;
        e2=e2-find_e;
        
        q_ins=sqrt((e1*q2^2-e2*q1^2)/(e1-e2));
        
        Q3=q_ins^2/(q1^2-q2^2);
        Q1=q1^2/(q1^2-q2^2);
        Q2=q2^2/(q1^2-q2^2);
        
        tv_e1=v_ep1(t_index-1);
        tv_e2=v_ep1(t_index);
        
        v_ins=Q3*(tv_e1-tv_e2)+Q1*tv_e2-Q2*tv_e1;
        
        length_vk1_ins=length_vk1_ins+1;
        v_k1_ins(length_vk1_ins)=q_ins;
        v_ep1_ins(length_vk1_ins)=v_ins;
        
    end
    
end

num_vk1=length(v_k1);

for i=1:length_vk1_ins
    
    tq=v_k1_ins(i);
    tv=v_ep1_ins(i);
    
    t_index=find(v_k1>tq,1);
    v_k1=[v_k1(1:t_index-1),tq,v_k1(t_index:num_vk1)];
    v_ep1=[v_ep1(1:t_index-1),tv,v_ep1(:,t_index:num_vk1)];
    num_vk1=num_vk1+1;
    
end

%%%%%%%%%%%%%%%%%%%%

for i_ins=1:NUM_int
    
    find_e=v_int_point(i_ins);  
    tvq=v_k2;
    tve=v_ep2;
    
    v_index=find((tve(1:end-1)-find_e).*(tve(2:end)-find_e)<0);
    
    for i2=1:length(v_index)
        
        t_index=v_index(i2)+1;
        
        q1=tvq(t_index-1);
        q2=tvq(t_index);
        
        e1=tve(t_index-1);
        e2=tve(t_index);
        
        e1=e1-find_e;
        e2=e2-find_e;
        
        q_ins=sqrt((e1*q2^2-e2*q1^2)/(e1-e2));
        
        Q3=q_ins^2/(q1^2-q2^2);
        Q1=q1^2/(q1^2-q2^2);
        Q2=q2^2/(q1^2-q2^2);
        
        tv_e1=v_ep2(t_index-1);
        tv_e2=v_ep2(t_index);
        
        v_ins=Q3*(tv_e1-tv_e2)+Q1*tv_e2-Q2*tv_e1;
        
        length_vk2_ins=length_vk2_ins+1;
        v_k2_ins(length_vk2_ins)=q_ins;
        v_ep2_ins(length_vk2_ins)=v_ins;
        
    end
    
end

num_vk2=length(v_k2);

for i=1:length_vk2_ins
    
    tq=v_k2_ins(i);
    tv=v_ep2_ins(i);
    
    t_index=find(v_k2>tq,1);
    v_k2=[v_k2(1:t_index-1),tq,v_k2(t_index:num_vk2)];
    v_ep2=[v_ep2(1:t_index-1),tv,v_ep2(:,t_index:num_vk2)];
    num_vk2=num_vk2+1;
    
end

Num_k1=num_vk1;
Num_k2=num_vk2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   we have just finished the inserting points, 
%   next we prepare for other thermodynamic quantities

v_length=v_k1(2:Num_k1)-v_k1(1:Num_k1-1);
v_k1_nt=zeros(1,Num_k1);
v_k1_nt(1:Num_k1-1)=1/2*v_length;
v_k1_nt(2:Num_k1)=v_k1_nt(2:Num_k1)+1/2*v_length;

v_length=v_k2(2:Num_k2)-v_k2(1:Num_k2-1);
v_k2_nt=zeros(1,Num_k2);
v_k2_nt(1:Num_k2-1)=1/2*v_length;
v_k2_nt(2:Num_k2)=v_k2_nt(2:Num_k2)+1/2*v_length;

[tmk,tmq]=meshgrid(v_k2,v_k2);
m_a2_22=a2(tmk-tmq)+a2(tmk+tmq);

[tmk,tmq]=meshgrid(v_k1,v_k2);
m_a1_21=a1(tmk-tmq)+a1(tmk+tmq);  %  2 transform to 1

[tmk,tmq]=meshgrid(v_k2,v_k1);
m_a1_12=a1(tmk-tmq)+a1(tmk+tmq); %  1 transfrom to 2

v_f_ep1=fn_F(v_ep1);
v_f_ep2=fn_F(v_ep2);

v_g_ep1=fn_G(v_ep1);
v_g_ep2=fn_G(v_ep2);

v_h_ep1=fn_H(v_ep1);
v_h_ep2=fn_H(v_ep2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    Let's start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  now: dmu

v_ep1_dmu_ori=-1*ones(1,Num_k1);
v_ep2_dmu_ori=-2*ones(1,Num_k2);

v_ep1_dmu=v_ep1_dmu_ori;
v_ep2_dmu=v_ep2_dmu_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dmu)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dmu_pre=v_ep1_dmu;
    v_ep2_dmu_pre=v_ep2_dmu;
    
    v_ep1_dmu=v_ep1_dmu_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dmu)*m_a1_21;
    v_ep2_dmu=v_ep2_dmu_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dmu)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dmu)*m_a2_22;
    
    flag1=max(abs(v_ep1_dmu-v_ep1_dmu_pre));
    flag2=max(abs(v_ep2_dmu-v_ep2_dmu_pre));
    flag=max(flag1,flag2);
    
end

result(2)=-1/pi/2*sum(v_k1_nt.*v_g_ep1.*v_ep1_dmu)...
    -1/pi*sum(v_k2_nt.*v_g_ep2.*v_ep2_dmu);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   next: dh

v_ep1_dh_ori=-1/2*ones(1,Num_k1);
v_ep2_dh_ori=0*ones(1,Num_k2);

v_ep1_dh=v_ep1_dh_ori;
v_ep2_dh=v_ep2_dh_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dh)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dh_pre=v_ep1_dh;
    v_ep2_dh_pre=v_ep2_dh;
    
    v_ep1_dh=v_ep1_dh_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dh)*m_a1_21;
    v_ep2_dh=v_ep2_dh_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dh)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dh)*m_a2_22;
    
    flag1=max(abs(v_ep1_dh-v_ep1_dh_pre));
    flag2=max(abs(v_ep2_dh-v_ep2_dh_pre));
    flag=max(flag1,flag2);
    
end

result(3)=-1/pi/2*sum(v_k1_nt.*v_g_ep1.*v_ep1_dh)...
    -1/pi*sum(v_k2_nt.*v_g_ep2.*v_ep2_dh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   next: dmudmu

v_ep1_dmudmu_ori=-(v_k2_nt.*v_h_ep2.*v_ep2_dmu.^2)*m_a1_21;
v_ep2_dmudmu_ori=-(v_k1_nt.*v_h_ep1.*v_ep1_dmu.^2)*m_a1_12...
                -(v_k2_nt.*v_h_ep2.*v_ep2_dmu.^2)*m_a2_22;

v_ep1_dmudmu=v_ep1_dmudmu_ori;
v_ep2_dmudmu=v_ep2_dmudmu_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dmudmu)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dmudmu_pre=v_ep1_dmudmu;
    v_ep2_dmudmu_pre=v_ep2_dmudmu;
    
    v_ep1_dmudmu=v_ep1_dmudmu_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dmudmu)*m_a1_21;
    v_ep2_dmudmu=v_ep2_dmudmu_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dmudmu)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dmudmu)*m_a2_22;
    
    flag1=max(abs(v_ep1_dmudmu-v_ep1_dmudmu_pre));
    flag2=max(abs(v_ep2_dmudmu-v_ep2_dmudmu_pre));
    flag=max(flag1,flag2);
    
end

result(4)=-1/pi/2*sum(v_k1_nt.*(v_g_ep1.*v_ep1_dmudmu+v_h_ep1.*v_ep1_dmu.^2))...
    -1/pi*sum(v_k2_nt.*(v_g_ep2.*v_ep2_dmudmu+v_h_ep2.*v_ep2_dmu.^2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  next: dhdh

v_ep1_dhdh_ori=-(v_k2_nt.*v_h_ep2.*v_ep2_dh.^2)*m_a1_21;
v_ep2_dhdh_ori=-(v_k1_nt.*v_h_ep1.*v_ep1_dh.^2)*m_a1_12...
                -(v_k2_nt.*v_h_ep2.*v_ep2_dh.^2)*m_a2_22;

v_ep1_dhdh=v_ep1_dhdh_ori;
v_ep2_dhdh=v_ep2_dhdh_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dhdh)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dhdh_pre=v_ep1_dhdh;
    v_ep2_dhdh_pre=v_ep2_dhdh;
    
    v_ep1_dhdh=v_ep1_dhdh_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dhdh)*m_a1_21;
    v_ep2_dhdh=v_ep2_dhdh_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dhdh)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dhdh)*m_a2_22;
    
    flag1=max(abs(v_ep1_dhdh-v_ep1_dhdh_pre));
    flag2=max(abs(v_ep2_dhdh-v_ep2_dhdh_pre));
    flag=max(flag1,flag2);
    
end

result(5)=-1/pi/2*sum(v_k1_nt.*(v_g_ep1.*v_ep1_dhdh+v_h_ep1.*v_ep1_dh.^2))...
    -1/pi*sum(v_k2_nt.*(v_g_ep2.*v_ep2_dhdh+v_h_ep2.*v_ep2_dh.^2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% prepare for dT
v_h1_ep1=1/T*(v_f_ep1-v_g_ep1.*v_ep1);
v_h1_ep2=1/T*(v_f_ep2-v_g_ep2.*v_ep2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  next: dT

v_ep1_dT_ori=-(v_k2_nt.*v_h1_ep2)*m_a1_21;
v_ep2_dT_ori=-(v_k1_nt.*v_h1_ep1)*m_a1_12...
                -(v_k2_nt.*v_h1_ep2)*m_a2_22;

v_ep1_dT=v_ep1_dT_ori;
v_ep2_dT=v_ep2_dT_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dT)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dT_pre=v_ep1_dT;
    v_ep2_dT_pre=v_ep2_dT;
    
    v_ep1_dT=v_ep1_dT_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dT)*m_a1_21;
    v_ep2_dT=v_ep2_dT_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dT)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dT)*m_a2_22;
    
    flag1=max(abs(v_ep1_dT-v_ep1_dT_pre));
    flag2=max(abs(v_ep2_dT-v_ep2_dT_pre));
    flag=max(flag1,flag2);
    
end

result(6)=-1/pi/2*sum(v_k1_nt.*(v_g_ep1.*v_ep1_dT+v_h1_ep1))...
    -1/pi*sum(v_k2_nt.*(v_g_ep2.*v_ep2_dT+v_h1_ep2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% prepare for dTdT
v_h2_ep1=v_h_ep1.*(v_ep1/T-v_ep1_dT).^2;
v_h2_ep2=v_h_ep2.*(v_ep2/T-v_ep2_dT).^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  next: dTdT

v_ep1_dTdT_ori=-(v_k2_nt.*v_h2_ep2)*m_a1_21;
v_ep2_dTdT_ori=-(v_k1_nt.*v_h2_ep1)*m_a1_12...
                -(v_k2_nt.*v_h2_ep2)*m_a2_22;

v_ep1_dTdT=v_ep1_dTdT_ori;
v_ep2_dTdT=v_ep2_dTdT_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dTdT)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dTdT_pre=v_ep1_dTdT;
    v_ep2_dTdT_pre=v_ep2_dTdT;
    
    v_ep1_dTdT=v_ep1_dTdT_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dTdT)*m_a1_21;
    v_ep2_dTdT=v_ep2_dTdT_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dTdT)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dTdT)*m_a2_22;
    
    flag1=max(abs(v_ep1_dTdT-v_ep1_dTdT_pre));
    flag2=max(abs(v_ep2_dTdT-v_ep2_dTdT_pre));
    flag=max(flag1,flag2);
    
end

result(7)=-1/pi/2*sum(v_k1_nt.*(v_g_ep1.*v_ep1_dTdT+v_h2_ep1))...
    -1/pi*sum(v_k2_nt.*(v_g_ep2.*v_ep2_dTdT+v_h2_ep2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% prepare for dTdmu
v_h3_ep1=-v_h_ep1.*v_ep1_dmu.*(v_ep1/T-v_ep1_dT);
v_h3_ep2=-v_h_ep2.*v_ep2_dmu.*(v_ep2/T-v_ep2_dT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  next: dTdmu

v_ep1_dTdmu_ori=-(v_k2_nt.*v_h3_ep2)*m_a1_21;
v_ep2_dTdmu_ori=-(v_k1_nt.*v_h3_ep1)*m_a1_12...
                -(v_k2_nt.*v_h3_ep2)*m_a2_22;

v_ep1_dTdmu=v_ep1_dTdmu_ori;
v_ep2_dTdmu=v_ep2_dTdmu_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dTdmu)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dTdmu_pre=v_ep1_dTdmu;
    v_ep2_dTdmu_pre=v_ep2_dTdmu;
    
    v_ep1_dTdmu=v_ep1_dTdmu_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dTdmu)*m_a1_21;
    v_ep2_dTdmu=v_ep2_dTdmu_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dTdmu)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dTdmu)*m_a2_22;
    
    flag1=max(abs(v_ep1_dTdmu-v_ep1_dTdmu_pre));
    flag2=max(abs(v_ep2_dTdmu-v_ep2_dTdmu_pre));
    flag=max(flag1,flag2);
    
end

result(8)=-1/pi/2*sum(v_k1_nt.*(v_g_ep1.*v_ep1_dTdmu+v_h3_ep1))...
    -1/pi*sum(v_k2_nt.*(v_g_ep2.*v_ep2_dTdmu+v_h3_ep2));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% prepare for dTdh
v_h4_ep1=-v_h_ep1.*v_ep1_dh.*(v_ep1/T-v_ep1_dT);
v_h4_ep2=-v_h_ep2.*v_ep2_dh.*(v_ep2/T-v_ep2_dT);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  next: dTdh

v_ep1_dTdh_ori=-(v_k2_nt.*v_h4_ep2)*m_a1_21;
v_ep2_dTdh_ori=-(v_k1_nt.*v_h4_ep1)*m_a1_12...
                -(v_k2_nt.*v_h4_ep2)*m_a2_22;

v_ep1_dTdh=v_ep1_dTdh_ori;
v_ep2_dTdh=v_ep2_dTdh_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dTdh)*m_a1_12;

flag_delta=0.1^10;
flag=1;
while flag>flag_delta

    v_ep1_dTdh_pre=v_ep1_dTdh;
    v_ep2_dTdh_pre=v_ep2_dTdh;
    
    v_ep1_dTdh=v_ep1_dTdh_ori-(v_k2_nt.*v_g_ep2.*v_ep2_dTdh)*m_a1_21;
    v_ep2_dTdh=v_ep2_dTdh_ori-(v_k1_nt.*v_g_ep1.*v_ep1_dTdh)*m_a1_12 ...
        -(v_k2_nt.*v_g_ep2.*v_ep2_dTdh)*m_a2_22;
    
    flag1=max(abs(v_ep1_dTdh-v_ep1_dTdh_pre));
    flag2=max(abs(v_ep2_dTdh-v_ep2_dTdh_pre));
    flag=max(flag1,flag2);
    
end

result(9)=-1/pi/2*sum(v_k1_nt.*(v_g_ep1.*v_ep1_dTdh+v_h4_ep1))...
    -1/pi*sum(v_k2_nt.*(v_g_ep2.*v_ep2_dTdh+v_h4_ep2));


% DO NOT FORGET TO TIMES 2
%%%%%%%%%%%%%%%%%%
result=result*2;%%
%%%%%%%%%%%%%%%%%%

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%    END of Main program    %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function f=fn_F(x)

global T

tt=-x./T;
ind1=find(tt>33);
ind2=find(tt<=33);

f(ind1)=x(ind1);
f(ind2)=-T*log(1+exp(tt(ind2)));

end

function f=fn_G(x)

global T
f=1./(1+exp(x/T));

end

function f=fn_H(x)

global T
f=-1./(exp(x/T)+exp(-x/T)+2)./T;

end