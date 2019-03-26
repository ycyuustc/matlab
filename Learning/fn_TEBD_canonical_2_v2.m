function [s1_new,t1_new,s2_new,t2_new,etaR,etaL] = fn_TEBD_canonical_2_v2(s1,t1,s2,t2)

TR = ncon({t1,s2,t2,s1,conj(t1),conj(s2),conj(t2),conj(s1)},...
    {[-1,1,7],[1,2],[2,3,8],[3,-3],...
    [-2,4,7],[4,5],[5,6,8],[6,-4]});
TL = ncon({s1,t1,s2,t2,conj(s1),conj(t1),conj(s2),conj(t2)},...
    {[-3,1],[1,2,7],[2,3],[3,-1,8],...
    [-4,4],[4,5,7],[5,6],[6,-2,8]});
ts2 = size(TR); tsR = [ts2,ones(1,4-length(ts2))];
mR = reshape(TR,tsR(1)*tsR(2),tsR(3)*tsR(4));
ts2 = size(TL); tsL = [ts2,ones(1,4-length(ts2))];
mL = reshape(TL,tsL(1)*tsL(2),tsL(3)*tsL(4));

[vR,etaR] = eigs(mR,1); vR = vR/sign(vR(1));
mX = reshape(vR,tsR(3),tsR(4));
[v,d] = eig(mX);
X = v*diag(sqrt(diag(d)))*v';
inv_X = v*diag(1./sqrt(diag(d)))*v';

[vL,etaL] = eigs(mL,1); vL = vL/sign(vL(1));
mY = reshape(vL,tsR(3),tsR(4));
[v,d] = eig(mY);
Y = v*diag(sqrt(diag(d)))*v';
inv_Y = v*diag(1./sqrt(diag(d)))*v';
Y = transpose(Y);
inv_Y = transpose(inv_Y);

tm = Y*s1*X;
[U,S,V] = svd(tm,'econ'); V = V';

t1_new = ncon({V*inv_X,t1},{[-1,1],[1,-2,-3]});
t2_new = ncon({t2,inv_Y*U},{[-1,1,-3],[1,-2]});
s1_new = S;
s2_new = s2;

end