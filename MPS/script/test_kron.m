sz=[1,0;0,-1];
sx=[0,1;1,0];
sy=[0,-1j;1j,0];

dt=0.05;

tic;
h=kron(sx,sx)+kron(sy,sy)+kron(sz,sz);
w=expm(-1j*dt*h);
w=reshape(w,[2,2,2,2]);
w=permute(w,[1,3,2,4]);
w=reshape(w,[4,4]);
[U,S,V]=svd2(w);
eta=size(S,1);
U=U*sqrt(S);
V=sqrt(S)*V;
U=reshape(U,[2,2,eta]);
U=permute(U,[4,3,2,1]);
V=reshape(V,[eta,2,2]);
V=permute(V,[1,4,3,2]);
toc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 上面是文章里的算法。我下面的算法和上面的算法等价。
tic;
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
toc;
%  至此 mu 和 U 相等， mv 和 v 相等。完成分解。
%  但上面的算法效率稍高。





