N=5;

Num=40;

vh=linspace(0.2,8,Num);
venergy=zeros(1,Num);
vmps=cell(1,Num);
v_1b=cell(1,Num);
m_2b=cell(1,Num);

for ih=1:Num
    
    h=vh(ih);
    [E0,mps]=fn_simpleheisenberg(N,h);
    venergy(ih)=E0;
    vmps{ih}=mps;
    [v_1b{ih},m_2b{ih}]=fn_datamining_mps(mps);
    
end

figure(4);
plot(vh,venergy,'*');

%%%%%%%%%%% 计算每个site上的极化 %%%%%%%%%%
sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1j;1j,0];
sz=sz/2;
sx=sx/2;
sy=sy/2;
id=eye(2);

mpola=zeros(Num,N);
for ih=1:Num
    t=v_1b{ih};
    for j=1:N
        onebody=t{j};
        pola=fn_contract(onebody,2,[1,2],sz,2,[1,2]);
        mpola(ih,j)=pola;
    end
    figure(9);
    hold on;
    plot(vh(ih)*ones(1,N),mpola(ih,:),'*');
end

figure(5);
plot(vh,sum(mpola,2),'*');

%%%%%%%%% 计算totalspin^2 %%%%%%%%%%

hspin=cell(3,2);
hspin{1,1}=sx;
hspin{1,2}=sx;
hspin{2,1}=sy;
hspin{2,2}=sy;
hspin{3,1}=sz;
hspin{3,2}=sz;

v_totalspin2=zeros(1,Num);

for ih=1:Num
    
    t=m_2b{ih};
    totalspin2=N*3/4;
    for i=1:N
        for j=1:N
            
            if i~=j
                for ispin=1:3
                    net=t{i,j};
                    net=fn_contract(net,4,[1,2],hspin{ispin,1},2,[1,2]);
                    net=fn_contract(net,2,[1,2],hspin{ispin,2},2,[1,2]);
                    totalspin2=totalspin2+net;
                end
            end
        end
    end
    
    v_totalspin2(ih)=totalspin2;
    
end

figure(6);
plot(vh,v_totalspin2,'*');


