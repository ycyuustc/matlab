

op = rand(2,2);
tm = fn_generate_mpo(op,3,2);

lambda = 0.01;
Num = 4;

Sz = [1/2,0;0,-1/2];
Sp = [0,1;0,0];
Sm = [0,0;1,0];
SI = eye(2);

cSz = cell(1,Num);
cSp = cell(1,Num);
cSm = cell(1,Num);
cSI = cell(1,Num);


for i=1:Num
   cSz{i} = fn_generate_mpo(Sz,Num,i);
   cSp{i} = fn_generate_mpo(Sp,Num,i);
   cSm{i} = fn_generate_mpo(Sm,Num,i);
end

hset = cell(1,Num);
gset = cell(1,Num);
for i=1:Num-1
    hset{i} = cSz{i}*cSz{i+1} + 0.5*(cSp{i}*cSm{i+1}+cSm{i}*cSp{i+1});
    gset{i} = expm(-lambda*hset{i});
end
hset{Num} = cSz{Num}*cSz{1} + 0.5*(cSp{Num}*cSm{1}+cSm{Num}*cSp{1});
gset{Num} = expm(-lambda*hset{Num});

H = 0;
for i = 1:Num-1
    H = H + hset{i};
end
H_cut = H - hset{4};

H_cut_L = 0;
for i=1:(Num/2-1)
    H_cut_L = H_cut_L + hset{i};
end
H_cut_L = reshape(H_cut_L,[2^(Num/2),2^(Num/2),2^(Num/2),2^(Num/2)]);
H_cut_L = ncon({H_cut_L},{[1,-1,1,-2]});

H_cut_R = 0;
for i=(Num/2):(Num-1)
    H_cut_R = H_cut_R + hset{i};
end

opts.disp  = 0; opts.issym = 1; opts.real  = 1;
[v_eig_cut,Energy_cut] = eigs(1/2*(H_cut+H_cut'),1,'SA',opts);
[v_eig,Energy] = eigs(1/2*(H+H'),1,'SA',opts);

psi_cut_L = reshape(v_eig_cut,2^(Num/2),2^(Num/2));
psi_cut_L = transpose(psi_cut_L);
rho_cut_L = psi_cut_L*psi_cut_L';

[V_cut_L,D_cut_L] = eig(rho_cut_L);
[D_cut_L, Index_cut_L] = sort(diag(D_cut_L),'descend');
V_cut_L = V_cut_L(:,Index_cut_L);

psi = reshape(v_eig,2^(Num/2),2^(Num/2));
psi = transpose(psi);
rho = psi*psi';

[V,D] = eig(rho);
[D, Index] = sort(diag(D),'descend');
V = V(:,Index);

tr = expm(-lambda*H);
% tm1 = gset{1}*gset{2}*gset{3}*gset{4}*gset{5}*gset{6};
% tm2 = gset{1}*gset{3}*gset{5}*gset{2}*gset{4}*gset{6};

%%%%%%%%%%%%%%%%%%%%%%%%\

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

cS1 = cell(1,3);
cS1{1} = Sz;
cS1{2} = Sp;
cS1{3} = Sm;

cS2 = cS1;

cBsL = cS1;
cBsR = cS1;
cBeL = cS1;
cBeR = cS1;

h = @(x,y) expm(-lambda*(x{1}*y{1} + 0.5*(x{2}*y{3} + x{3}*y{2})));

cBsR_super = cell(1,3);
cBeR_super = cell(1,3);
cS1_super = cell(1,3);
cS2_super = cell(1,3);
for i = 1:10
    
    BsI = eye(NKeep);
    BeI = eye(NKeep);
    
    cS1{1} = kron(Sz,BsI);
    cS1{2} = kron(Sp,BsI);
    cS1{3} = kron(Sm,BsI);
    cS1_super{1} = kron(kron(kron(Sz,BsI),SI),BeI);
    cS1_super{2} = kron(kron(kron(Sp,BsI),SI),BeI);
    cS1_super{3} = kron(kron(kron(Sm,BsI),SI),BeI);
    
    
    cBsL{1} = kron(SI,cBsL{1});
    cBsL{2} = kron(SI,cBsL{2});
    cBsL{3} = kron(SI,cBsL{3}); 
    Ts = kron(SI,Ts);
    cBsR{1} = kron(SI,cBsR{1});
    cBsR{2} = kron(SI,cBsR{2});
    cBsR{3} = kron(SI,cBsR{3});    
    Ts = h(cS1,cBsL)*Ts;
    
    cS2{1} = kron(Sz,BeI);
    cS2{2} = kron(Sp,BeI);
    cS2{3} = kron(Sm,BeI);
    cS2_super{1} = kron(kron(kron(SI,BsI),Sz),BeI);
    cS2_super{2} = kron(kron(kron(SI,BsI),Sp),BeI);
    cS2_super{3} = kron(kron(kron(SI,BsI),Sm),BeI);
    
    cBeL{1} = kron(SI,cBeL{1});
    cBeL{2} = kron(SI,cBeL{2});
    cBeL{3} = kron(SI,cBeL{3}); 
    Te = kron(SI,Te);
    cBeR{1} = kron(SI,cBeR{1});
    cBeR{2} = kron(SI,cBeR{2});
    cBeR{3} = kron(SI,cBeR{3});    
    Te = h(cS2,cBeL)*Te;
    
    Ts_super = kron(kron(Ts,SI),BeI);
    Te_super = kron(kron(SI,BsI),Te);
    
    cBsR_super{1} = kron(kron(cBsR{1},SI),BeI);
    cBsR_super{2} = kron(kron(cBsR{2},SI),BeI);
    cBsR_super{3} = kron(kron(cBsR{3},SI),BeI);
    
    cBeR_super{1} = kron(kron(SI,BsI),cBeR{1});
    cBeR_super{2} = kron(kron(SI,BsI),cBeR{2});
    cBeR_super{3} = kron(kron(SI,BsI),cBeR{3});
    
    T = Ts_super*h(cBsR_super,cS2_super)*Te_super*h(cBeR_super,cS1_super);
     
end


function res = fn_generate_mpo(op,Num,site)

Num_left = site-1;
Num_right = Num-site;

res = kron(kron(eye(2^Num_left),op),eye(2^Num_right));

end


