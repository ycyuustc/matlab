Sz = [1/2,0;0,-1/2];
Sp = [0,1;0,0];
Sm = [0,0;1,0];
SI = eye(2);
dim = 2;
m = 20;
m_id2 = eye(2);

f = @(x,y) reshape(ncon({x,y},{[-2,-4],[-1,-3]}),...
    size(x,1)*size(y,1),size(x,2)*size(y,2));
f4 = @(x1,x2,x3,x4) f(f(f(x1,x2),x3),x4);

Hmid = f(Sz,Sz)+0.5*(f(Sp,Sm)+f(Sm,Sp));



BH = zeros(dim,dim);
% BI = eye(dim);
Bz = Sz;
Bp = Sp;
Bm = Sm;

BHr = zeros(dim,dim);
% BIr = eye(dim);
Bzr = Sz;
Bpr = Sp;
Bmr = Sm;

NKeep = dim;

fprintf('Iter\tTotalLength\t Energy per site \t\t\t\t Trunc\n');

for i = 1:100
   
   % prepare the operators: 
   BI = eye(NKeep);
   BIr = eye(NKeep);
      
   p2 = f4(BI,Sp,SI,BIr);
   m2 = f4(BI,Sm,SI,BIr);
   z2 = f4(BI,Sz,SI,BIr);
   
   p3 = f4(BI,SI,Sp,BIr);
   m3 = f4(BI,SI,Sm,BIr);
   z3 = f4(BI,SI,Sz,BIr);
   
   BH = f(BH,SI) + f(Bz,Sz) + 0.5*(f(Bp,Sm) + f(Bm,Sp));
   Bz = f(BI,Sz); Bp = f(BI,Sp); Bm = f(BI,Sm);
   BH_super = f(f(BH,SI),BIr);
   BHr = f(SI,BHr) + f(Sz,Bzr) + 0.5*(f(Sm,Bpr) + f(Sp,Bmr));
   Bzr = f(Sz,BI); Bpr = f(Sp,BI); Bmr = f(Sm,BI);
   BHr_super = f(f(BI,SI),BHr);
   
   HSB = BH_super...
       + z2*z3 + 0.5*(p2*m3 + m2*p3)...
       + BHr_super;
   
   HSB = (HSB+HSB')/2;
   
   opts.disp = 0; opts.issym = 1; opts.real = 1;
   [Psi, Energy] = eigs(HSB,1,'SA',opts);
   
   len = length(Psi);
   PsiMatrix = reshape(Psi,sqrt(len),sqrt(len));
   PsiMatrix = transpose(PsiMatrix);
   
   Rho = PsiMatrix*PsiMatrix';
   [V,D] = eig(Rho);
   [D,Index] = sort(diag(D),'descend');
   V = V(:,Index);
   NKeep = min(size(D,1),m);
   T = V(:,1:NKeep);
   
   Rhor = PsiMatrix'*PsiMatrix;
   [Vr,Dr] = eig(Rhor);
   [Dr,Indexr] = sort(diag(Dr),'descend');
   Vr = Vr(:,Indexr);
   NKeepr = min(size(Dr,1),m);
   Tr = Vr(:,1:NKeepr);
   
   TruncationError = 1-sum(D(1:NKeep));
   
   BH = T'*BH*T;
   Bz = T'*Bz*T; 
   Bp = T'*Bp*T; 
   Bm = T'*Bm*T;
   
   BHr = Tr'*BHr*Tr;
   Bzr = Tr'*Bzr*Tr;
   Bpr = Tr'*Bpr*Tr; 
   Bmr = Tr'*Bmr*Tr;
   
   fprintf('%d\t\t %d\t\t\t %f\t\t\t %f\n',i, 2*i+2, Energy/2/(i+1), TruncationError);
end