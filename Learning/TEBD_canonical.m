T = rand(4,4,4);
s = diag(rand(1,4));

T2 = ncon({T,conj(T)},{[-1,-3,1],[-2,-4,1]});

TL = ncon({s,conj(s),T2},{[-1,1],[-2,2],[1,2,-3,-4]});
TR = ncon({T2,s,conj(s)},{[-1,-2,1,2],[1,-3],[2,-4]});

ts = size(TR);
mR = reshape(TR,ts(1)*ts(2),ts(3)*ts(4));
[v,~] = eigs(mR,1); v = v/sign(v(1));
tm = reshape(v,ts(1),ts(2));
[v,d] = eig(tm);
X = v*sqrt(d)*v';
inv_X = v*diag(1./sqrt(diag(d)))*v';
disp(ncon({TR,X,conj(X)},{[-1,-2,1,2],[1,3],[2,3]})./(X*X'));

TL = permute(TL,[3,4,1,2]);
ts = size(TL);
mL = reshape(TL,ts(1)*ts(2),ts(3)*ts(4));
[v,eta] = eigs(mL,1); v = v/sign(v(1));
tm = reshape(v,ts(1),ts(2));
[v,d] = eig(tm);
Y = v*sqrt(d)*v';
inv_Y = v*diag(1./sqrt(diag(d)))*v';
Y = transpose(Y);
inv_Y = transpose(inv_Y);
TL = permute(TL,[3,4,1,2]);
disp(ncon({Y,conj(Y),TL},{[1,3],[2,3],[1,2,-1,-2]})./(transpose(Y)*transpose(Y)'));

tm = Y*s*X;
[u,s_new,v] = svd(tm,'econ'); v = v';
T_new = ncon({v*inv_X,T,inv_Y*u},{[-1,1],[1,2,-3],[2,-2]});

disp('check');
disp(ncon({T,s},{[1,2,-1],[2,1]})-ncon({T_new,s_new},{[1,2,-1],[2,1]}));
T2_new = ncon({T_new,conj(T_new)},{[-1,-3,1],[-2,-4,1]});
TL_new = ncon({s_new,conj(s_new),T2_new},{[-1,1],[-2,2],[1,2,-3,-4]});
TR_new = ncon({T2_new,s_new,conj(s_new)},{[-1,-2,1,2],[1,-3],[2,-4]});

disp(ncon({TR_new},{[-1,-2,1,1]}));
disp(ncon({TL_new},{[1,1,-1,-2]}));

% T_new, s_new finished!

