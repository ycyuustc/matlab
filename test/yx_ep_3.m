M = 20;
Num = 2*M;

a=1.2;
b=0.4;
c=1;

tv1 = reshape([a*ones(1,M);c*ones(1,M)],[1,Num]);
tv2 = reshape([b*ones(1,M);c*ones(1,M)],[1,Num]);

m = diag(tv1(1:end-1),-1) + diag(tv2(1:end-1),1);
% m(1,Num) = c; m(Num,1) = c;  m(1,1) = V;
disp('H='); disp(num2str(m));
[v,d] = eig(m);