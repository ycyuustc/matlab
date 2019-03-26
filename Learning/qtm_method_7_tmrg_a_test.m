beta = 1;
M = 80;
lambda = beta/M*2;
m = 16;

figure(2);
hold on;


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

local_horizontal = 2*(s1z*s2z + 1/2*(s1p*s2m+s1m*s2p))-0.5*eye(4);
tm = expm(-lambda*local_horizontal);
a = 1/2*(tm(1,1)+tm(4,4));
b = 1/2*(tm(2,2)+tm(3,3));
c = 1/2*(tm(2,3)+tm(3,2));
cS1 = cell(1,4);
cS1{1} = kron(Sz,SI);
cS1{2} = kron(Sp,SI);
cS1{3} = kron(Sm,SI);
cS1{4} = kron(SI,SI);
cS2 = cell(1,4);
cS2{1} = kron(SI,Sz);
cS2{2} = kron(SI,Sp);
cS2{3} = kron(SI,Sm);
cS2{4} = kron(SI,SI);

g = @(x,y) (c+b)/2*x{4}*y{4}+2*(c-b)*x{1}*y{1}+a*(x{2}*y{3}+x{3}*y{2});
T12 = reshape(g(cS1,cS2),[2,2,2,2]);
T23 = T12;
T34 = T23;
T41 = T34;
T_ori = T12;

dim = 2;
NKeep = dim;

T_super = mcon({T12,T23,T34,T41},...
    {[3,4,-7,-8],[-2,-3,2,3],[1,2,-5,-6],[-4,-1,4,1]});
T_super = 0.5*(T_super+permute(T_super,[7,6,5,8,3,2,1,4]));
m_super = reshape(T_super,[NKeep*dim*NKeep*dim,NKeep*dim*NKeep*dim]);
[psi_R,Energy_R] = eigs(m_super,1);
psi_R = reshape(psi_R,[NKeep,dim,NKeep,dim]);
psi_L = permute(psi_R,[3,2,1,4]);
norm_factor = sqrt(mcon({psi_R,psi_L},{[1,2,3,4],[1,2,3,4]}));
psi_R = psi_R/norm_factor;
psi_L = psi_L/norm_factor;
rho_s = mcon({psi_R,psi_L},{[2,-1,-2,1],[2,-3,-4,1]});
rho_s = reshape(rho_s,[NKeep*dim,NKeep*dim]);
NKeep_pre = NKeep;
NKeep = min(NKeep*dim,m);
[V_R,~] = eigs(rho_s,NKeep);
[V_L,~] = eigs(transpose(rho_s),NKeep);
[P,D] = eig(transpose(V_L)*V_R);
D = diag(1./sqrt(diag(D)));
VR = V_R*P*D;
VL = V_L*transpose(inv(P))*D;
U = reshape(VR,dim,NKeep_pre,NKeep);
V = reshape(VL,dim,NKeep_pre,NKeep);

for i=1:(M/2-2)
    if mod(i,2) == 1
        T_super = mcon({T12,T23,T_ori,T_ori,T34,T41},...
            {[5,6,-11,-12],[-4,-5,4,5],[3,4,-9,-10],...
            [-2,-3,2,3],[1,2,-7,-8],[-6,-1,6,1]});
        T_super = mcon({T_super,U,V,U,V},...
            {[1,2,-2,3,4,-4,5,6,-6,7,8,-8],[2,1,-1],[3,4,-3],...
            [7,8,-7],[6,5,-5]});
        T23 = mcon({T23,T_ori,U,V},...
            {[1,2,3,4],[-1,3,-3,5],[5,4,-4],[1,2,-2]});
        T34 = mcon({T34,T_ori,U,V},...
            {[1,3,4,5],[2,-2,3,-4],[2,1,-1],[5,4,-3]});
    else
        T_super = mcon({T12,T23,T_ori,T_ori,T34,T41},...
            {[5,6,-11,-12],[4,-5,-10,5],[-3,-4,3,4],...
            [2,3,-8,-9],[1,-2,-7,2],[-6,-1,6,1]});
        T_super = mcon({T_super,U,V,U,V},...
            {[1,2,-2,3,4,-4,5,6,-6,7,8,-8],[2,1,-1],[3,4,-3],...
            [7,8,-7],[6,5,-5]});
        T23 = mcon({T23,T_ori,U,V},...
            {[3,2,4,5],[-1,1,-3,3],[4,5,-4],[1,2,-2]});
        T34 = mcon({T34,T_ori,U,V},...
            {[1,2,3,4],[4,-2,5,-4],[2,1,-1],[5,3,-3]});
    end
    T12 = mcon({T12,eye(dim),U,V},...
        {[2,-2,4,-4],[1,3],[3,4,-3],[1,2,-1]});
    T41 = mcon({T41,eye(dim),U,V},...
        {[-1,1,-3,3],[2,4],[2,1,-2],[4,3,-4]});
    
    T_super = 0.5*(T_super+permute(T_super,[7,6,5,8,3,2,1,4]));  
    m_super = reshape(T_super,[NKeep*dim*NKeep*dim,NKeep*dim*NKeep*dim]);
    
    [psi_R,Energy_R] = eigs(m_super,1);
    disp(sqrt(Energy_R));
    plot(2*i+4,real(Energy_R),'bo');
    
    psi_R = reshape(psi_R,[NKeep,dim,NKeep,dim]);
    psi_L = permute(psi_R,[3,2,1,4]);
    norm_factor = sqrt(mcon({psi_R,psi_L},{[1,2,3,4],[1,2,3,4]}));
    %    disp(norm_factor);
    psi_R = psi_R/norm_factor;
    psi_L = psi_L/norm_factor;
    
    rho_s = mcon({psi_R,psi_L},{[2,-1,-2,1],[2,-3,-4,1]});
    rho_s = reshape(rho_s,[NKeep*dim,NKeep*dim]);
    
    NKeep_pre = NKeep;
    NKeep = min(NKeep*dim,m);
%     [V_R,~] = eigs(rho_s,NKeep);
%     [V_L,~] = eigs(transpose(rho_s),NKeep);
%     [P,D] = eig(transpose(V_L)*V_R);
%     D = diag(1./sqrt(diag(D)));
%     VR = V_R*P*D;
%     VL = V_L*transpose(inv(P))*D;

[V_R,D_R] = eig(rho_s);
[~,Index] = sort(abs(diag(D_R)), 'descend');
V_R = V_R(:,Index);
V_L = inv(V_R);
T_R = V_R(:,1:NKeep);
T_L = V_L(1:NKeep,:);
VL = transpose(T_L);
VR = T_R;

%     [V_R,D0,V_L] = eig(rho_s);
%     [~,Index] = sort(abs(diag(D0)), 'descend');
% V_R = V_R(:,Index);
% V_L = V_L(:,Index);
%     [P,D] = eig(transpose(V_L)*V_R);
%     D = diag(1./sqrt(diag(D)));
%     VR = V_R*P*D;
%     VL = V_L*transpose(inv(P))*D;


% [V_R,~] = eigs(rho_s,NKeep);
% [V_L,~] = eigs(transpose(rho_s),NKeep);
% [P,D] = eig(transpose(V_L)*V_R);
% D = diag(1./sqrt(diag(D)));
% VR = V_R*P*D;
% VL = V_L*transpose(inv(P))*D;
    
    U = reshape(VR,dim,NKeep_pre,NKeep);
    V = reshape(VL,dim,NKeep_pre,NKeep);
    
    T12 = 0.5*(T12 + permute(T41,[4,3,2,1]));
    T41 = 0.5*(T41 + permute(T12,[4,3,2,1]));
    T23 = 0.5*(T23 + permute(T34,[4,3,2,1]));
    T34 = 0.5*(T34 + permute(T23,[4,3,2,1]));

end