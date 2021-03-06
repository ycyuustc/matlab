beta = 0.1;
M = 20;
lambda = beta/M*2;
m = 16;

Sz = [1/2 0;0 -1/2];    %sigma_z
Sp = [0 1;0 0];         %s^+�� Sp=Sx+Im*Sy;   
Sm = [0 0;1 0];         %s^-�� Sm=Sx-Im*Sy;  
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

dim = 2;
g = @(x,y) (c+b)/2*x{4}*y{4}+2*(c-b)*x{1}*y{1}+a*(x{2}*y{3}+x{3}*y{2});
TSL = reshape(g(cS1,cS2),[2,2,2,2]);    
TSR = TSL;
TEL = TSL;
TER = TSL;
TS = eye(dim);
TE = eye(dim);

NKeep = dim;
for i=1:(M/2-1)
    
   if mod(i,2) == 1      
        T_super = ncon({TSL,TS,TSR,TEL,TE,TER},...
            {[5,6,-7,-8],[4,5],[-2,-3,3,4],[2,3,-5,-6],[1,2],[-4,-1,6,1]});
   else   
       TSLR = 0.5*ncon({TSL,TSR},{[1,-3,-5,-6],[-1,-2,-4,1]});
       TSLR = TSLR + 0.5*ncon({TSL,TSR},{[-2,-3,1,-6],[-1,1,-4,-5]});
       TELR = 0.5*ncon({TEL,TER},{[1,-3,-5,-6],[-1,-2,-4,1]});
       TELR = TELR + 0.5*ncon({TEL,TER},{[-2,-3,1,-6],[-1,1,-4,-5]});
       T_super = ncon({TSLR,TS,TE,TELR},...
           {[3,4,-7,-8],[-3,3],[1,-5],[-4,-1,1,4]});
         
   end
    
   m_super = reshape(T_super,[NKeep*dim*NKeep*dim,NKeep*dim*NKeep*dim]);
   
   [psi_R,Energy_R] = eigs(m_super,1);
   disp(sqrt(Energy_R));
   
   psi_R = reshape(psi_R,[NKeep,dim,NKeep,dim]);
   psi_L = permute(psi_R,[3,2,1,4]);
   norm_factor = sqrt(ncon({psi_R,psi_L},{[1,2,3,4],[1,2,3,4]}));
%    disp(norm_factor);
   psi_R = psi_R/norm_factor;
   psi_L = psi_L/norm_factor;
   
   rho_s = ncon({psi_R,psi_L},{[2,-1,-2,1],[2,-3,-4,1]});
   rho_s = reshape(rho_s,[NKeep*dim,NKeep*dim]);
   
   NKeep_pre = NKeep;
   NKeep = min(NKeep*dim,m);
   [V_R,D_R] = eig(rho_s);
   [~,Index] = sort(abs(diag(D_R)), 'descend');
   V_R = V_R(:,Index);
   V_L = inv(V_R);
   T_R = V_R(:,1:NKeep);
   T_L = V_L(1:NKeep,:);
   VL = transpose(T_L);
   VR = T_R;
   
   U = reshape(VR,dim,NKeep_pre,NKeep);
   V = reshape(VL,dim,NKeep_pre,NKeep);
   
   T12 = ncon({T12,eye(dim),U,V},...
       {[2,-2,4,-4],[1,3],[3,4,-3],[1,2,-1]});   
   T41 = ncon({T41,eye(dim),U,V},...
       {[-1,1,-3,3],[2,4],[2,1,-2],[4,3,-4]});
   
   if mod(i,2)==1
       T23 = ncon({T23,T_ori,U,V},...
           {[1,2,3,4],[-1,3,-3,5],[5,4,-4],[1,2,-2]});    
       T34 = ncon({T34,T_ori,U,V},...
           {[1,3,4,5],[2,-2,3,-4],[2,1,-1],[5,4,-3]}); 
   else
       T23 = ncon({T23,T_ori,U,V},...
           {[3,2,4,5],[-1,1,-3,3],[4,5,-4],[1,2,-2]});     
       T34 = ncon({T34,T_ori,U,V},...
           {[1,2,3,4],[4,-2,5,-4],[2,1,-1],[5,3,-3]});       
   end
     
end




