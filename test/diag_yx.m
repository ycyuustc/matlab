b = 0.8;
c = 1/b;

Num = 40;

vb = b*ones(1,Num-1);
vc = c*ones(1,Num-1);
va = 2*rand(1,Num);

mh = diag(va) + diag(vb,1) + diag(vc,-1);
mh(Num,1) = b;
mh(1,Num) = c;

[v,d] = eig(mh);
disp(max(max(abs(v*(d/v)-mh))));