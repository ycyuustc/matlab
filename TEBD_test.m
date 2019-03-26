sx = [0,1;1,0];
sI = eye(2);
sz = [1,0;0,-1];
hxx = kron(sx,sx);
thxx = reshape(hxx,[2,2,2,2]);
hzz = kron(sz,sz);
hz = kron(sz,sI) + kron(sI,sz);
hx = kron(sx,sI) + kron(sI,sx);
hz1 = kron(sz,sI);
hz2 = kron(sI,sz);
thz1 = reshape(hz1,[2,2,2,2]);
thz2 = reshape(hz2,[2,2,2,2]);
thzz = reshape(hzz,[2,2,2,2]);

epsilon = -0.01;
% mh = -1/2*(hxx + 1.01/2*hz);
mh = -hzz + 0.01*hz;
h12 = reshape(mh,[2,2,2,2]);
hm = expm(epsilon*mh);
hamil = reshape(hm,[2,2,2,2]);

TA = rand(2,2,2);
TB = rand(2,2,2);
sAB = diag(rand(1,2));
sBA = diag(rand(1,2));

D = 64;
delta = 0.01;
num_walk = 1000;

[TA,sAB,TB,sBA] = fn_iTEBD_ori(TA,sAB,TB,sBA,h12,D,delta,num_walk);
