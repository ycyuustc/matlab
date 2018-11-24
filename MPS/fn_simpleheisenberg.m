function [E0,mps]=fn_simpleheisenberg(N,h)

precision=1e-5;
D=8;

M=3*(N-1);
% M=M+N;
hset=cell(M,N);
hset_h=cell(N,N);

sx=[0,1;1,0];
sy=[0,-1j;1j,0];
sz=[1,0;0,-1];

id=eye(2);

for m=1:M
    for j=1:N
        hset{m,j}=id;
    end
end

for m=1:N
    for j=1:N
        hset_h{m,j}=id;
    end
end

for j=1:N
   
    hset_h{j,j}=-h*sz;
    
end

for j=1:(N-1)
   
    hset{3*(j-1)+1,j}=sx;
    hset{3*(j-1)+1,j+1}=sx;
    hset{3*(j-1)+2,j}=sy;
    hset{3*(j-1)+2,j+1}=sy;
    hset{3*(j-1)+3,j}=sz;
    hset{3*(j-1)+3,j+1}=sz;

end

hset=[hset;hset_h];

[E0,mps]=fn_minimizeE(hset,D,precision);

end