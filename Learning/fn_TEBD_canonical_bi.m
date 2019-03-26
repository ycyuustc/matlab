function [A_new,lambda_new,B_new,theta_new,sR,sL,etaR,etaL] = fn_TEBD_canonical_bi(A,lambda,B,theta)

TR = ncon({A,lambda,conj(B),conj(theta)},{[-1,2,1],[2,-3],[-2,3,1],[3,-4]});
ts2=size(TR); ts=[ts2,ones(1,4-length(ts2))];
mR = reshape(TR,ts(1)*ts(2),ts(3)*ts(4));
[v,etaR] = eigs(mR,1); v = v/sign(v(1));
tm = reshape(v,ts(1),ts(2));
[u,s,v] = svd(tm,'econ');
X = u;
U = v;
sR = s;

TL = ncon({lambda,A,conj(theta),conj(B)},{[-3,2],[2,-1,1],[-4,3],[3,-2,1]});
ts2=size(TL); ts=[ts2,ones(1,4-length(ts2))];
mL = reshape(TL,ts(1)*ts(2),ts(3)*ts(4));
[v,etaL] = eigs(mL,1); v = v/sign(v(1));
tm = reshape(v,ts(1),ts(2));
[u,s,v] = svd(tm,'econ');
Y = conj(u);
V = conj(v);
sL = s;

A_new = ncon({X',A,Y},{[-1,1],[1,2,-3],[2,-2]});
B_new = ncon({U',B,V},{[-1,1],[1,2,-3],[2,-2]});
lambda_new = Y'*lambda*X;
theta_new = V'*theta*U;

end