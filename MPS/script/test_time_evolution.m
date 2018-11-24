N=10;
D=5;
precision=1e-5;
dt=0.03;
jflipped=5;

oset=cell(1,N);
for j=1:N
    oset{1,j}=id;
end
oset{1,jflipped}=sz;

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

mps0=cell(1,N);
for j=1:N
    if j==jflipped
        state=[0;1];
    else
        state=[1;0];
    end
    mps0{j}=reshape(state,[1,1,2]);
end

mps=mps0;

mzvalues=zeros(1,50);
for step=1:50
    fprintf('step %2d:',step);
    [mps,K1]=fn_reduceD(mps,mpo_even,D,precision);
    [mps,K2]=fn_reduceD(mps,mpo_odd,D,precision);
    mz=fn_expectationvalue(mps,oset);
    mzvalues(step)=mz;
   fprintf('mz=%g\n',mz);
end
