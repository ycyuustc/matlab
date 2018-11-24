N=20;
D=32;

Num=4;

precision=1e-5;

vh=linspace(0,4,Num);

venergy=zeros(1,100);
cmps=cell(1,100);

for ih=1:Num

h=vh(ih);
M=3*(N-1);
M=M+N;
hset=cell(M,N);

sx=[0,1;1,0];
sy=[0,-1j;1j,0];
sz=[1,0;0,-1];

id=eye(2);

for m=1:M
    for j=1:N
        hset{m,j}=id;
    end
end

for j=1:(N-1)
   
    hset{3*(j-1)+1,j}=sx;
    hset{3*(j-1)+1,j+1}=sx;
    hset{3*(j-1)+2,j}=sy;
    hset{3*(j-1)+2,j+1}=sy;
    hset{3*(j-1)+3,j}=sz;
    hset{3*(j-1)+3,j+1}=sz;

end

for j=1:N
    
    hset{3*N-3+j,j}=-h*sz;
    
end

[E0,mps0]=fn_minimizeE(hset,D,precision);

venergy(ih)=E0;
cmps{ih}=mps0;

disp(E0);

end