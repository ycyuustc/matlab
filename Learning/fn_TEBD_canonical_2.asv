function [s1_new,t1_new,s2_new,t2_new] = fn_TEBD_canonical_2(s1,t1,s2,t2)

s = s1;
T = ncon({t1,s2,t2},{[-1,1,-3],[1,2],[2,-2,-4]});
ts = size(T);
T = reshape(T,[ts(1),ts(2),ts(3)*ts(4)]);
[T,s] = fn_TEBD_canonical(T,s);
T = reshape(T,ts);
T = permute(T,[1,3,2,4]);
mT = reshape(T,ts(1)*ts(3),ts(2)*ts(4));
[U,S,V] = svd(mT,'econ'); V = V';
s1_new = s;
s2_new = S;
t1_new = reshape(U,[ts(1),ts(3),ts(2)*ts(4)]);
t1_new = permute(t1_new,[1,3,2]);
t2_new = reshape(V,[ts(1)*ts(3),ts(2),ts(4)]);

end