Sz = [1/2 0;0 -1/2];    %sigma_z
Sp = [0 1;0 0];         %s^+£¬ Sp=Sx+Im*Sy;   
Sm = [0 0;1 0];         %s^-£¬ Sm=Sx-Im*Sy;  
SI = eye(2);

s1z = kron(Sz,SI);
s2z = kron(SI,Sz);
s1p = kron(Sp,SI);
s2p = kron(SI,Sp);
s1m = kron(Sm,SI);
s2m = kron(SI,Sm);

local_1 = 2*(s1z*s2z + 1/2*(s1p*s2m+s1m*s2p))-1/2*eye(4);
disp(local_1);
local_2 = 1/2*s1z*s2z-1/4*(s1p*s2p+s1m*s2m)+1/2*(s1p*s2m+s1m*s2p)+1/8*eye(4);
disp(local_2);

beta = 0.01;
Num = 4;
M = 4;
lambda = beta/M;
m = 12;

local_exp_1 = expm(-lambda*local_1);
tm = local_exp_1;
a = 1/2*(tm(1,1)+tm(4,4));
b = 1/2*(tm(2,2)+tm(3,3));
c = 1/2*(tm(2,3)+tm(3,2));


cSz = cell(1,Num);
cSp = cell(1,Num);
cSm = cell(1,Num);
cSI = cell(1,Num);

for i=1:Num
   cSz{i} = fn_generate_mpo(Sz,Num,i);
   cSp{i} = fn_generate_mpo(Sp,Num,i);
   cSm{i} = fn_generate_mpo(Sm,Num,i);
   cSI{i} = fn_generate_mpo(SI,Num,i);
end

hset = cell(1,Num);
gset = cell(1,Num);
for i=1:Num-1
    hset{i} = cSz{i}*cSz{i+1} + 0.5*(cSp{i}*cSm{i+1}+cSm{i}*cSp{i+1});
    gset{i} = (a+b)/2*cSI{i}*cSI{i+1} + 2*(a-b)*cSz{i}*cSz{i+1} ...
        + c*(cSp{i}*cSm{i+1} + cSm{i}*cSp{i+1});
end
hset{Num} = cSz{Num}*cSz{1} + 0.5*(cSp{Num}*cSm{1}+cSm{Num}*cSp{1});
gset{Num} = (a+b)/2*cSI{Num}*cSI{1} + 2*(a-b)*cSz{Num}*cSz{1} ...
        + c*(cSp{Num}*cSm{1} + cSm{Num}*cSp{1});

transverse_matrix_1 = 1;
transverse_matrix_2 = 1;
for i=1:2:(Num-1)
    transverse_matrix_1 = transverse_matrix_1*gset{i};
    transverse_matrix_2 = transverse_matrix_2*gset{i+1};
end

Z1 = trace((transverse_matrix_1*transverse_matrix_2)^(M/2));
disp(['Z1=', num2str(Z1)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cSz = cell(1,M);
cSp = cell(1,M);
cSm = cell(1,M);
cSI = cell(1,M);

for i=1:M
   cSz{i} = fn_generate_mpo(Sz,M,i);
   cSp{i} = fn_generate_mpo(Sp,M,i);
   cSm{i} = fn_generate_mpo(Sm,M,i);
   cSI{i} = fn_generate_mpo(SI,M,i);
end

hset2 = cell(1,M);
gset2 = cell(1,M);
for i=1:M-1
    hset2{i} = 0.5*cSz{i}*cSz{i+1} + 0.5*(cSp{i}*cSm{i+1}+cSm{i}*cSp{i+1})...
        -0.25*(cSp{i}*cSp{i+1}+cSm{i}*cSm{i+1}) + 0.125*cSI{i}*cSI{i+1};
    gset2{i} = a/2*cSI{i}*cSI{i+1} + 2*a*cSz{i}*cSz{i+1} ...
        + b*(cSp{i}*cSp{i+1} + cSm{i}*cSm{i+1}) ...
        + c*(cSp{i}*cSm{i+1} + cSm{i}*cSp{i+1});
end
hset2{M} = 0.5*cSz{M}*cSz{1} + 0.5*(cSp{M}*cSm{1}+cSm{M}*cSp{1})...
    -0.25*(cSp{M}*cSp{1} + cSm{M}*cSm{1}) + 0.125*cSI{M}*cSI{1};
gset2{M} = a/2*cSI{M}*cSI{1} + 2*a*cSz{M}*cSz{1} ...
        + b*(cSp{M}*cSp{1} + cSm{M}*cSm{1}) ...
        + c*(cSp{M}*cSm{1} + cSm{M}*cSp{1});
horizontal_matrix_1 = 1;
horizontal_matrix_2 = 1;
for i=1:2:(M-1)
    horizontal_matrix_1 = horizontal_matrix_1*gset2{i};
    horizontal_matrix_2 = horizontal_matrix_2*gset2{i+1};
end

Z2 = trace((horizontal_matrix_1*horizontal_matrix_2)^(Num/2));
disp(['Z2=', num2str(Z2)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cSz = cell(1,M);
cSp = cell(1,M);
cSm = cell(1,M);
cSI = cell(1,M);

for i=1:M
   cSz{i} = fn_generate_mpo(Sz,M,i);
   cSp{i} = fn_generate_mpo(Sp,M,i);
   cSm{i} = fn_generate_mpo(Sm,M,i);
   cSI{i} = fn_generate_mpo(SI,M,i);
end

hset3 = cell(1,M);
gset3 = cell(1,M);
for i=1:M-1
    gset3{i} = (c+b)/2*cSI{i}*cSI{i+1} + 2*(c-b)*cSz{i}*cSz{i+1} ...
        + a*(cSp{i}*cSm{i+1} + cSm{i}*cSp{i+1});
end
gset3{M} = (c+b)/2*cSI{M}*cSI{1} + 2*(c-b)*cSz{M}*cSz{1} ...
        + a*(cSp{M}*cSm{1} + cSm{M}*cSp{1});

horizontal_matrix_1 = 1;
horizontal_matrix_2 = 1;
for i=1:2:(M-1)
    horizontal_matrix_1 = horizontal_matrix_1*gset3{i};
    horizontal_matrix_2 = horizontal_matrix_2*gset3{i+1};
end

Z3 = trace((horizontal_matrix_1*horizontal_matrix_2)^(Num/2));
disp(['Z3=', num2str(Z3)]);
disp(eig(horizontal_matrix_1*horizontal_matrix_2));
tT = horizontal_matrix_2*horizontal_matrix_1;
tT2 = horizontal_matrix_1*horizontal_matrix_2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tm = (c+b)/2*eye(4) + 2*(c-b)*s1z*s2z + a*(s1p*s2m+s1m*s2p);
aa = tm2(1,1); bb = tm2(2,2); cc = tm2(2,3);

cS1 = cell(1,4);
cS1{1} = Sz;
cS1{2} = Sp;
cS1{3} = Sm;
cS1{4} = SI;

cS2 = cS1;

cBsL = cS1;
cBsR = cS1;
cBeL = cS1;
cBeR = cS1;

cS1_super = cell(1,4);
cS2_super = cell(1,4);
cBsL_super = cell(1,4);
cBsR_super = cell(1,4);
cBeL_super = cell(1,4);
cBeR_super = cell(1,4);
cBsL_block = cell(1,4);
cBsR_block = cell(1,4);
cBeL_block = cell(1,4);
cBeR_block = cell(1,4);
cS2s_block = cell(1,4);
cS2e_block = cell(1,4);

g = @(x,y) (c+b)/2*x{4}*y{4} + 2*(c-b)*x{1}*y{1} ...
        + a*(x{2}*y{3} + x{3}*y{2});
glog = @(x,y) (aa+bb)/2*x{4}*y{4} + 2*(aa-bb)*x{1}*y{1} ...
        + cc*(x{2}*y{3} + x{3}*y{2});
f4 = @(x1,x2,x3,x4) kron(kron(kron(x1,x2),x3),x4);
f2 = @(x1,x2) kron(x1,x2);
    
dim = 2;
NKeep = dim;
BsLz = Sz;
BsLp = Sp;
BsLm = Sm;
Ts = eye(dim);

BeLz = Sz;
BeLp = Sp;
BeLm = Sm;
Te = eye(dim);

for i=1:20
    
    BsI = eye(NKeep);
    BeI = eye(NKeep);
    disp(eigs(Ts));
    disp(eigs(Te));
    
    % construction of the super block and block
    
    for k = 1:4      
        cS1_super{k} = f4(cS1{k},BsI,SI,BeI);
        
        cBsL_super{k} = f4(SI,cBsL{k},SI,BeI);
        cBsR_super{k} = f4(SI,cBsR{k},SI,BeI);
           
        cS2_super{k} = f4(SI,BsI,cS2{k},BeI);
        
        cBeL_super{k} = f4(SI,BsI,SI,cBeL{k});
        cBeR_super{k} = f4(SI,BsI,SI,cBeR{k});
        
        cBsL_block{k} = f2(cBsL{k},SI);
        cBsR_block{k} = f2(cBsR{k},SI);
        
        cS2s_block{k} = f2(BsI,cS2{k});
            
        cBeL_block{k} = f2(SI,cBeL{k});
        cBeR_block{k} = f2(SI,cBeR{k});
        
        cS2e_block{k} = f2(cS2{k},BeI);
          
    end
    cS1_super{4} = eye(size(cS1_super{4}));
    cBsL_super{4} = eye(size(cBsL_super{4}));
    cBsR_super{4} = eye(size(cBsR_super{4}));
    cBeL_super{4} = eye(size(cBeL_super{4}));
    cBeR_super{4} = eye(size(cBeR_super{4}));
    
    cBsL_block{4} = eye(size(cBsL_block{4}));
    cBsR_block{4} = eye(size(cBsR_block{4}));
    cS2s_block{4} = eye(size(cS2s_block{4}));
    cBeL_block{4} = eye(size(cBeL_block{4}));
    cBeR_block{4} = eye(size(cBeR_block{4}));
    cS2e_block{4} = eye(size(cS2e_block{4}));
    
    Ts_super = f4(SI,Ts,SI,BeI); 
    Te_super = f4(SI,BsI,SI,Te);  
    Ts_block = f2(Ts,SI);
    Te_block = f2(SI,Te);
       
    g12 = g(cS1_super,cBsL_super);
    g12_log = glog(cS1_super,cBsL_super);
    g23 = g(cBsR_super,cS2_super);
    g23_log = glog(cBsR_super,cS2_super);
    g34 = g(cS2_super,cBeL_super);
    g34_log = glog(cS2_super,cBeL_super);
    g41 = g(cBeR_super,cS1_super);
    g41_log = glog(cBeR_super,cS1_super);
    
    g23_block = g(cBsR_block,cS2s_block);
    g34_block = g(cS2e_block,cBeL_block);
     
    if mod(i,2) == 0
         T_super2 = g41*g34*Ts_super*Te_super*g23*g12;
         T_super3 = g34*g41*Ts_super*Te_super*g23*g12;
        T_super = expm(g41_log+g34_log)*Ts_super*Te_super...
            *expm(g23_log+g12_log);
        Ts_block = Ts_block*g23_block;
        Te_block = g34_block*Te_block;
    else
         T_super2 = g41*g23*Ts_super*Te_super*g34*g12;
         T_super3 = g23*g41*Ts_super*Te_super*g34*g12;
        T_super = expm(g41_log+g23_log)*Ts_super*Te_super...
            *expm(g34_log+g12_log);
        Ts_block = g23_block*Ts_block;
        Te_block = Te_block*g34_block;
    end
    
    disp(max(max(abs(T_super2-T_super))));
    
%     disp(eig(T_super));
    
%     T_super = reshape(T_super,[NKeep,dim,NKeep,dim,NKeep,dim,NKeep,dim]);
    
    %% %%%%%%%% transition of the operators %%%%%%%%%%%%%%%%%%%%%

    [psi_R,Energy_R] = eigs(T_super,1);
%     [psi_L,Energy_L] = eigs(transpose(T_super),1);
    
%     norm_factor = sqrt(sum((psi_L).*psi_R));
%     disp(norm_factor);
%     psi_R = psi_R/norm_factor;
%     psi_L = psi_L/norm_factor;
%     
    disp(Energy_R);
%     disp(Energy_L);
    
    psi_R = reshape(psi_R,[NKeep,dim,NKeep,dim]);
%     psi_L = reshape(psi_L,[NKeep,dim,NKeep,dim]);   
    psi_L = permute(psi_R,[3,2,1,4]);
%     psi_R = 1/2*(psi_R + permute(psi_R,[3,4,1,2]));
%     psi_L = 1/2*(psi_L + permute(psi_L,[3,4,1,2]));
%      psi_R2 = permute(psi_R,[1,4,3,2]);
%      disp(psi_R2-psi_L);
%     psi_L = psi_R2;
%     
    norm_factor = sqrt(ncon({psi_R,psi_L},{[1,2,3,4],[1,2,3,4]}));
    disp(norm_factor);
    psi_R = psi_R/norm_factor;
    psi_L = psi_L/norm_factor;
    
    rho_s = ncon({psi_R,psi_L},{[2,-1,-2,1],[2,-3,-4,1]});
    rho_s = reshape(rho_s,[NKeep*dim,NKeep*dim]);
    rho_e = ncon({psi_R,psi_L},{[-1,-2,1,2],[-3,-4,1,2]});
    rho_e = reshape(rho_e,[NKeep*dim,NKeep*dim]);
    
    disp(trace(rho_s));
    disp(trace(rho_e));
    
    NKeep_pre = NKeep;
    NKeep = min(NKeep*dim,m);
    
    %% %%%%%%%% the system block operators %%%%%%%%%%%%%
    [V_R,D_R] = eig(rho_s);
    [~,Index] = sort(abs(diag(D_R)), 'descend');
    V_R = V_R(:,Index);
    V_L = inv(V_R);
    T_R = V_R(:,1:NKeep);
    T_L = V_L(1:NKeep,:);

%     [V_R,D,V_L] = eig(rho_s);
%     [~,Index] = sort(abs(diag(D)),'descend');
%     V_R = V_R(:,Index);
%     V_L = V_L(:,Index);
%     
%     T_L = V_L(:,1:NKeep);
%     T_R = V_R(:,1:NKeep);
%     
    
%     [V_L,D_L] = eigs(transpose(rho_s),NKeep);
%     [~,Index] = sort(abs(diag(D_L)), 'descend');
%     V_L = V_L(:,Index);
%     T_L = V_L(:,1:NKeep);

    
    
%     [u,s,v] = svd(transpose(T_L)*T_R);
%     VR = T_R*v*diag(1./sqrt(diag(s)));
%     VL = T_L*conj(u)*diag(1./sqrt(diag(s)));
%     tm = transpose(T_L)*T_R;
%     tm = diag(tm);
%     tm = sqrt(tm);
%     tm = 1./tm;
%     tm = diag(tm);
%     VR = T_R*tm;
%     VL = T_L*tm;
     VL = transpose(T_L);
     VR = T_R;
    
    Ts = transpose(VL)*Ts_block*VR;
    for k=1:4
       cBsL{k} = transpose(VL)*cBsL_block{k}*VR;
       cBsR{k} = transpose(VL)*cS2s_block{k}*VR; 
    end
    cBsL{4} = eye(size(cBsL{4}));
    cBsR{4} = eye(size(cBsR{4}));
    
    %% %%%%%%%% the enviroment block operators %%%%%%%%%
    
    tVR = reshape(VR,[dim,NKeep_pre,NKeep]);
    tVR = permute(tVR,[2,1,3]);
    tVR = reshape(tVR,dim*NKeep_pre,NKeep);
    
    tVL = reshape(VL,[dim,NKeep_pre,NKeep]);
    tVL = permute(tVL,[2,1,3]);
    tVL = reshape(tVL,dim*NKeep_pre,NKeep);
    
    VL = tVR;
    VR = tVL;
    
%     [V_R,D_R] = eig(rho_e);
%     [D_R,Index] = sort(abs(diag(D_R)), 'descend');
%     V_R = V_R(:,Index);
%     T_L = inv(V_R);
%     T_R = V_R(:,1:NKeep);
%     T_L = T_L(1:NKeep,:);
    
%     [V_R,D,V_L] = eig(rho_e);
%     [~,Index] = sort(abs(diag(D)),'descend');
%     V_R = V_R(:,Index);
%     V_L = V_L(:,Index);
%     
%     T_L = V_L(:,1:NKeep);
%     T_R = V_R(:,1:NKeep);
%     
%     [V_L,D_L] = eig(transpose(rho_e));
%     [D_L,Index] = sort(abs(diag(D_L)), 'descend');
%     V_L = V_L(:,Index);
%     T_L = V_L(:,1:NKeep);
%     [T_R,~] = eigs(rho_e,NKeep);
%     [T_L,~] = eigs(transpose(rho_e),NKeep);
%     
%     [u,s,v] = svd(transpose(T_L)*T_R);
%     VR = T_R*v*diag(1./sqrt(diag(s)));
%     VL = T_L*conj(u)*diag(1./sqrt(diag(s)));

%     tm = transpose(T_L)*T_R;
%     tm = diag(tm);
%     tm = sqrt(tm);
%     tm = 1./tm;
%     tm = diag(tm);
%     VR = T_R*tm;
%     VL = T_L*tm;
%      VR = T_R;
%      VL = transpose(T_L);
    
    Te = transpose(VL)*Te_block*VR;
    for k=1:4
       cBeL{k} = transpose(VL)*cS2e_block{k}*VR;
       cBeR{k} = transpose(VL)*cBeR_block{k}*VR; 
    end
    cBeL{4} = eye(size(cBeL{4}));
    cBeR{4} = eye(size(cBeR{4}));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    Ts = 1/2*(Ts + transpose(Te));
    Te = transpose(Ts);
    cBsL{1} = 1/2*(cBsL{1} + transpose(cBeR{1}));
    cBsL{2} = 1/2*(cBsL{2} + transpose(cBeR{3}));
    cBsL{3} = 1/2*(cBsL{3} + transpose(cBeR{2}));
    cBsR{1} = 1/2*(cBsR{1} + transpose(cBeL{1}));
    cBsR{2} = 1/2*(cBsR{2} + transpose(cBeL{3}));
    cBsR{3} = 1/2*(cBsR{3} + transpose(cBeL{2}));
    
    cBeR{1} = transpose(cBsL{1});
    cBeR{2} = transpose(cBsL{3});
    cBeR{3} = transpose(cBsL{2});
    cBeL{1} = transpose(cBsR{1});
    cBeL{2} = transpose(cBsR{3});
    cBeL{3} = transpose(cBsR{2});

    disp(max(max(abs(Ts - transpose(Te)))));
    disp(max(max(abs(cBsL{1} - transpose(cBeR{1})))));
    disp(max(max(abs(cBsL{2} - transpose(cBeR{3})))));
    disp(max(max(abs(cBsL{3} - transpose(cBeR{2})))));
    disp(max(max(abs(cBsR{1} - transpose(cBeL{1})))));
    disp(max(max(abs(cBsR{2} - transpose(cBeL{3})))));
    disp(max(max(abs(cBsR{3} - transpose(cBeL{2})))));

    disp(fn_judge(cBsL{1},cBsR{1}));
    
    disp(i);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function res = fn_judge(x,y)
    res = max(max(abs(x*y-y*x)));
end


function res = fn_generate_mpo(op,Num,site)

Num_left = site-1;
Num_right = Num-site;

res = kron(kron(eye(2^Num_left),op),eye(2^Num_right));

end



