sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1j;1j,0];
id=eye(2);

N=5;
D=5;
d=2;

% mpst=fn_createrandommps(N,D,d);
hset_sz=cell(N,N);
for i=1:N
    for j=1:N
        hset_sz{i,j}=reshape(id,[1,1,2,2]);
    end
end

for j=1:N
    hset_sz{j,j}=reshape(sz,[1,1,2,2]);
end

[e,n]=fn_expectationvalue(mpst,hset_sz);
% e是mz， n是norm。
% 以上是文献中的算法。下面我用datamining函数，能够得到一样的结果。
% 



[v_1b,m_2b]=fn_datamining_mps(mpst);

mz=0;
for j=1:N
    mj=v_1b{j};
    mzj=fn_contract(mj,2,[1,2],sz,2,[1,2]);
    mz=mz+mzj;
end

mz=mz/fn_contractmps(mpst);