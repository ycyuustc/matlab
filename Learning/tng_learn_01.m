N = 4;
M = 5;
D = 6;

B = rand(N,N,D)+1j*rand(N,N,D);
temp = ncon({B,conj(B)},{[1,2,3],[1,2,3]});
disp(temp);
B_norm = temp;
vt = rand(M,N,N);
temp = ncon({B,vt,conj(vt),conj(B)},{[1,2,6],[3,1,2],[3,4,5],[4,5,6]});
disp(temp);

temp = ncon({B,vt},{[1,2,-2],[-1,1,2]});
Bvnorm = ncon({temp,conj(temp)},{[1,2],[1,2]});
disp(Bvnorm);

Bd = conj(B);
vtd = conj(vt);

for i=1:10000
ev = ncon({B,conj(vt),conj(B)},{[-2,-3,3],[-1,1,2],[1,2,3]});
ev = reshape(ev,M,N*N);
[u,s,v] = svd(ev,'econ');
vt = v*u';
vt = reshape(vt,N,N,M);
vt = permute(vt,[3,1,2]);

temp = ncon({B,vt,conj(vt),conj(B)},{[1,2,6],[3,1,2],[3,4,5],[4,5,6]});
disp(B_norm-temp);
end




