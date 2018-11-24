N=3;
D=8;
precision=1e-4;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sx=1/sqrt(2)*[0,1,0;1,0,1;0,1,0];
sy=1j/sqrt(2)*[0,-1,0;1,0,-1;0,1,0];
sz=[1,0,0;0,0,0;0,0,-1];
sxx=sx*sx;
syy=sy*sy;
szz=sz*sz;
sxy=sx*sy;
syx=sy*sx;
sxz=sx*sz;
szx=sz*sx;
syz=sy*sz;
szy=sz*sy;

id=eye(3);

hset_h=cell(N,N);
hset_1=cell(3*(N-1),N);
hset_2=cell(9*(N-1),N);
hset_id=cell(1,N);

for j=1:N
   hset_id{1,j}=id; 
end


for i=1:N
    for j=1:N
        hset_h{j,i}=id;
    end
end

for i=1:N
    for j=1:3*(N-1)
        hset_1{j,i}=id;
    end
end

for i=1:N
    for j=1:9*(N-1)
        hset_2{j,i}=id;
    end
end

% for j=1:N
%    
%     hset_h{j,j}=sz;
%     
% end
% 
% for j=1:N-1
%    
%     hset_1{3*(j-1)+1,j}=sx;
%     hset_1{3*(j-1)+1,j+1}=sx;
%     
%     hset_1{3*(j-1)+2,j}=sy;
%     hset_1{3*(j-1)+2,j+1}=sy;
%     
%     hset_1{3*(j-1)+3,j}=sz;
%     hset_1{3*(j-1)+3,j+1}=sz;
%     
% end
% 
% 
% for j=1:(N-1)
%    
%     hset_2{9*(j-1)+1,j}=sxx;
%     hset_2{9*(j-1)+1,j+1}=sxx;
%     
%     hset_2{9*(j-1)+2,j}=syy;
%     hset_2{9*(j-1)+2,j+1}=syy;
%     
%     hset_2{9*(j-1)+3,j}=szz;
%     hset_2{9*(j-1)+3,j+1}=szz;
%     
%     hset_2{9*(j-1)+4,j}=sxy;
%     hset_2{9*(j-1)+4,j+1}=sxy;
%     
%     hset_2{9*(j-1)+5,j}=syx;
%     hset_2{9*(j-1)+5,j+1}=syx;
%     
%     hset_2{9*(j-1)+6,j}=sxz;
%     hset_2{9*(j-1)+6,j+1}=sxz;
%     
%     hset_2{9*(j-1)+7,j}=szx;
%     hset_2{9*(j-1)+7,j+1}=szx;
%     
%     hset_2{9*(j-1)+8,j}=syz;
%     hset_2{9*(j-1)+8,j+1}=syz;
%     
%     hset_2{9*(j-1)+9,j}=szy;
%     hset_2{9*(j-1)+9,j+1}=szy;
%    
%     
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% b2=1;
% b1=1;
h=0;

vbeta=linspace(0,2.0,10);
vy=zeros(1,10);
vy2=zeros(1,10);
vy3=zeros(1,10);
vy4=zeros(1,10);

for i=1:10

beta=vbeta(i);    
% b2=1/3*(1+1/2*x);
% b1=1/2*x;
% b0=1/3*(x-1);

b1=1;
b2=beta;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
b1=sqrt(b1);
b2=sqrt(b2);

% b1=-1;
% b2=-1;


for j=1:N
    
    hset_h{j,j}=sz*(-h);
    
end

for j=1:N-1
    
    hset_1{3*(j-1)+1,j}=-sx*b1;
    hset_1{3*(j-1)+1,j+1}=sx*b1;
    
    hset_1{3*(j-1)+2,j}=-sy*b1;
    hset_1{3*(j-1)+2,j+1}=sy*b1;
    
    hset_1{3*(j-1)+3,j}=-sz*b1;
    hset_1{3*(j-1)+3,j+1}=sz*b1;
    
end


for j=1:(N-1)
    
    hset_2{9*(j-1)+1,j}=-sxx*b2;
    hset_2{9*(j-1)+1,j+1}=sxx*b2;
    
    hset_2{9*(j-1)+2,j}=-syy*b2;
    hset_2{9*(j-1)+2,j+1}=syy*b2;
    
    hset_2{9*(j-1)+3,j}=-szz*b2;
    hset_2{9*(j-1)+3,j+1}=szz*b2;
    
    hset_2{9*(j-1)+4,j}=-sxy*b2;
    hset_2{9*(j-1)+4,j+1}=sxy*b2;
    
    hset_2{9*(j-1)+5,j}=-syx*b2;
    hset_2{9*(j-1)+5,j+1}=syx*b2;
    
    hset_2{9*(j-1)+6,j}=-sxz*b2;
    hset_2{9*(j-1)+6,j+1}=sxz*b2;
    
    hset_2{9*(j-1)+7,j}=-szx*b2;
    hset_2{9*(j-1)+7,j+1}=szx*b2;
    
    hset_2{9*(j-1)+8,j}=-syz*b2;
    hset_2{9*(j-1)+8,j+1}=syz*b2;
    
    hset_2{9*(j-1)+9,j}=-szy*b2;
    hset_2{9*(j-1)+9,j+1}=szy*b2;
    
    
end


hset=[hset_1;hset_2];


[E0,mps0]=fn_minimizeE(hset,D,precision);
disp(['E0=',num2str(E0),'   b2=',num2str(b2),'   b1=',num2str(b1)]);

vy(i)=(E0-2*(1-beta));
vy2(i)=E0;

vy3(i)=-0.5-2*beta-sqrt((2*beta-1)^2+5/4);
% vy4(i)=-sqrt((3*beta-1)^2+2^2);

end

