%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  此脚本旨在visualize那个spinon这种奇葩。
%  step 1: 解出一个反铁磁基态 mpsAFM
%  step 2: 翻转一个自旋，生成一个新的态 mpsEx
%  step 3: 计算这个态随时间的烟花，啦啦啦
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% step 1: generate ground state of mpsAFM
precision=1e-5;
D=6;
N=12;
h=0;
jflipped=5;
sz=[1,0;0,-1]/2;

[E0,mpsAFM]=fn_simpleheisenberg(N,h);

%%%%%%%%%% step 2: generate the excitation state 
sfan=[1,1;0,0];
id=eye(2);
sfan=reshape(sfan,[1,1,2,2]);
id=reshape(id,[1,1,2,2]);
mpoX=cell(1,N);
jflippped=5;
for j=1:N
    if j==jflipped
        mpoX{j}=sfan;
    else
        mpoX{j}=id;
    end
end

mpsEx=fn_reduceD(mpsAFM,mpoX,D,precision);

[v_1b,m_2b]=fn_datamining_mps(mpsAFM);
mpola0=zeros(1,N);
for j=1:N
    onebody=v_1b{j};
    pola=fn_contract(onebody,2,[1,2],sz,2,[1,2]);
    mpola0(j)=pola;
end

[v_1b,m_2b]=fn_datamining_mps(mpsEx);
mpolaex=zeros(1,N);
for j=1:N
    onebody=v_1b{j};
    pola=fn_contract(onebody,2,[1,2],sz,2,[1,2]);
    mpolaex(j)=pola;
end

%%%%%%%%%%% step 3: time evolution 
dt=0.01;

oset=cell(1,N);
for j=1:N
    oset{1,j}=id;
end
oset{1,jflipped}=sz;

% 重写，这里的sx,sy,sz实际上是\sigma_x,\sigma_y,\sigma_z
sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1j;1j,0];
id=eye(2);

sxx=fn_contract(sx,3,3,sx,3,3);
syy=fn_contract(sy,3,3,sy,3,3);
szz=fn_contract(sz,3,3,sz,3,3);

st=sxx+syy+szz;
st=permute(st,[1,3,2,4]);
st=reshape(st,[4,4]);

st=expm(-1j*dt*st);

st=reshape(st,[2,2,2,2]);
st=permute(st,[1,3,2,4]);
st=reshape(st,[4,4]);

[ms,mv,md]=svd2(st);
eta=size(mv,1);
ms=ms*sqrt(mv);
md=sqrt(mv)*md;

mu=reshape(ms,[2,2,eta]);
mv=reshape(md,[eta,2,2]);

mu=permute(mu,[4,3,2,1]);
mv=permute(mv,[1,4,3,2]);
mI=reshape(id,[1,1,2,2]);

mpo_even=cell(1,N);
mpo_odd=cell(1,N);
for j=1:N
    mpo_even{j}=mI;
    mpo_odd{j}=mI;
end

for j=1:2:(N-1)
    mpo_odd{j}=mu;
    mpo_odd{j+1}=mv;
end

for j=2:2:(N-1)
    mpo_even{j}=mu;
    mpo_even{j+1}=mv;
end

Num_step=1000;

mps=mpsEx; % 初态是这个spinon激发
msz=zeros(Num_step+1,N);

[v_1b,m_2b]=fn_datamining_mps(mps);
mpola=zeros(1,N);
for j=1:N
    onebody=v_1b{j};
    pola=fn_contract(onebody,2,[1,2],sz/2,2,[1,2]);
    mpola(j)=pola;
end
msz(1,:)=mpola;

for step=1:Num_step
    fprintf('step %2d:\n',step);
    [mps,K1]=fn_reduceD(mps,mpo_even,D,precision);
    [mps,K2]=fn_reduceD(mps,mpo_odd,D,precision);
    [v_1b,m_2b]=fn_datamining_mps(mps);
    mpola=zeros(1,N);
    for j=1:N
        onebody=v_1b{j};
        pola=fn_contract(onebody,2,[1,2],sz/2,2,[1,2]);
        mpola(j)=pola;
    end
    msz(step+1,:)=mpola;
    
end

msz=real(msz);

close(figure(12));
figure(12);
axis([-0.1,8.1,-1.2,1.2]);
for step=1:Num_step+1
    
    plot(msz(step,:),'*');
    axis([-0.1,N+1.1,-0.7,0.7]);
    pause(0.04);
end
   



