beta = 0.1;
M = 12;
h = 0;

m2id=eye(2);
m2z=[1,0;0,-1];
dT = 2;

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;

m4p2 = permute(m4p,[1,4,3,2]);

precision = 1e-8;

vD = 16;
vN = M/2;
vbeta = beta;

m_r = zeros(length(vD),length(vN),length(vbeta));
for ibeta = 1:length(vbeta)
    beta = vbeta(ibeta);    
    for iN = 1:length(vN)
        N = vN(iN);   
        for iD = 1:length(vD)            
            D = vD(iD);
            % ================================================
            hset = fn_generate_mpo(N*2,beta,h);
            % ================================================
             hset{1,1} = -hset{1,1};
%             hset{2,1} = hset{2,1};
            [E0,mps0]=fn_minimizeE(hset,D,precision);
            E0 = -E0;
            m_r(iD,iN,ibeta) = real(E0);  
            fprintf('max eigen value: %f \n',E0);
        end
    end
end



function res = fn_generate_mpo(M,beta,h)

N = M/2;
Z = [exp(beta*h/2),0;0,exp(-beta*h/2)];
dT = 2;

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;

m4p2 = permute(m4p,[1,4,3,2]);
lambda = beta/2/N;

p1 = 1/(1-lambda);
p2 = lambda/(lambda-1);

T1 = p1*m4p + p2*m4id;
T2 = p1*m4p2 + p2*m4id;

hset=cell(1,N);

bT1_left = fn_contract(Z,2,2,T1,4,1);
bT1_left = fn_contract(bT1_left,4,1,T1,4,1);
bT1_left = permute(bT1_left,[2,5,1,4,3,6]);
bT1_left = reshape(bT1_left,[1,dT^2,dT^2,dT^2]);

bT2_right = fn_contract(T2,4,3,T2,4,3);
bT2_right = permute(bT2_right,[1,4,2,5,3,6]);
bT2_right = reshape(bT2_right,[dT^2,1,dT^2,dT^2]);

hset{1,1} = bT1_left;
hset{1,N} = bT2_right;

for i=2:(N-1)
    if mod(i,2) == 1
        bT1 = fn_contract(T1,5,5,...
            reshape(T1,[1,dT,dT,dT,dT]),5,1);
        bT1 = permute(bT1,[1,5,3,7,2,6,4,8]);
        bT1 = reshape(bT1,[dT^2,dT^2,dT^2,dT^2]);
        hset{1,i} = bT1;
    else
        bT2 = fn_contract(T2,5,5,...
            reshape(T2,[1,dT,dT,dT,dT]),5,1);
        bT2 = permute(bT2,[1,5,3,7,2,6,4,8]);
        bT2 = reshape(bT2,[dT^2,dT^2,dT^2,dT^2]);
        hset{1,i} = bT2;
    end   
end

res = hset;
end
