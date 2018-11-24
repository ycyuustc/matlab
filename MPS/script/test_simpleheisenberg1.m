N=4;   % site 个数

Num=40; % 离散点个数

vh=linspace(0.2,5,Num);  % 外场取值范围
venergy=zeros(1,Num);
venergy1=zeros(1,Num);
vmps=cell(1,Num);
vmps1=cell(1,Num);
v_1b=cell(1,Num);
v_1b1=cell(1,Num);
m_2b=cell(1,Num);
m_2b1=cell(1,Num);

for ih=1:Num
    
    h=vh(ih);
    [E0,mps,E1,mps1]=fn_simpleheisenberg1(N,h);
    venergy(ih)=E0;
    venergy1(ih)=E1;
    vmps{ih}=mps;
    vmps1{ih}=mps1;
    [v_1b{ih},m_2b{ih}]=fn_datamining_mps(mps);
    [v_1b1{ih},m_2b1{ih}]=fn_datamining_mps(mps1);
    
end

close(figure(4));
figure(4);
plot(vh,venergy,'*');
ylabel('energy');
hold on;
plot(vh,venergy1,'r*');
ylabel('energy');


%%%%%%%%%%% 计算每个site上的极化 %%%%%%%%%%
sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1j;1j,0];
sz=sz/2;
sx=sx/2;
sy=sy/2;
id=eye(2);

close(figure(9));
close(figure(10));

mpola=zeros(Num,N);
mpola1=zeros(Num,N);
for ih=1:Num
    t=v_1b{ih};
    t1=v_1b1{ih};
    for j=1:N
        onebody=t{j};
        onebody1=t1{j};
        pola=fn_contract(onebody,2,[1,2],sz,2,[1,2]);
        pola1=fn_contract(onebody1,2,[1,2],sz,2,[1,2]);
        mpola(ih,j)=pola;
        mpola1(ih,j)=pola1;
    end
    
    figure(9);
    hold on;
    plot(vh(ih)*ones(1,N),mpola(ih,:),'*');
    ylabel('mz per site');
    
    figure(10);
    hold on;
    plot(vh(ih)*ones(1,N),mpola1(ih,:),'r*');
    ylabel('mz per site');
end
close(figure(5));
figure(5);
plot(vh,sum(mpola,2),'*');
ylabel('mz');
hold on;
plot(vh,sum(mpola1,2),'r*');
ylabel('mz');

%%%%%%%%% 计算totalspin^2 %%%%%%%%%%

hspin=cell(3,2);
hspin{1,1}=sx;
hspin{1,2}=sx;
hspin{2,1}=sy;
hspin{2,2}=sy;
hspin{3,1}=sz;
hspin{3,2}=sz;

v_totalspin2=zeros(1,Num);
v_totalspin21=zeros(1,Num);

for ih=1:Num
    
    t=m_2b{ih};
    t1=m_2b1{ih};
    totalspin2=N*3/4;
    totalspin21=N*3/4;
    for i=1:N
        for j=1:N
            
            if i~=j
                for ispin=1:3
                    net=t{i,j};
                    net1=t1{i,j};
                    net=fn_contract(net,4,[1,2],hspin{ispin,1},2,[1,2]);
                    net=fn_contract(net,2,[1,2],hspin{ispin,2},2,[1,2]);
                    net1=fn_contract(net1,4,[1,2],hspin{ispin,1},2,[1,2]);
                    net1=fn_contract(net1,2,[1,2],hspin{ispin,2},2,[1,2]);
                    totalspin2=totalspin2+net;
                    totalspin21=totalspin21+net1;
                end
            end
        end
    end
    
    v_totalspin2(ih)=totalspin2;
    v_totalspin21(ih)=totalspin21;
    
end

close(figure(6));
figure(6);
plot(vh,v_totalspin2,'*');
ylabel('s^2');
close(figure(7));
figure(7);
plot(vh,v_totalspin21,'r*');
ylabel('s^2');


