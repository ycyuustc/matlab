Sz = [1/2,0;0,-1/2];
Sp = [0,1;0,0];
Sm = [0,0;1,0];
dim = 2;
m = 50;
m_id2 = eye(2);

f = @(x,y) reshape(ncon({x,y},{[-2,-4],[-1,-3]}),...
    size(x,1)*size(y,1),size(x,2)*size(y,2));

Hmid = f(Sz,Sz)+0.5*(f(Sp,Sm)+f(Sm,Sp));

BlockH = zeros(dim,dim);
BlockSz = Sz;
BlockSp = Sp;
BlockSm = Sm;

BlockHr = zeros(dim,dim);
BlockSzr = Sz;
BlockSpr = Sp;
BlockSmr = Sm;

NKeep = dim;

for i = 1:100
    
   BlockI = eye(NKeep);
   
   BlockH = f(BlockH,m_id2) + f(BlockSz,Sz)...
       + 0.5*(f(BlockSp,Sm)+f(BlockSm,Sp));
   
   BlockSz = f(BlockI,Sz);
   BlockSp = f(BlockI,Sp);
   BlockSm = f(BlockI,Sm);
   
   BlockHr = f(m_id2,BlockHr) + f(Sz,BlockSzr)...
       +0.5*(f(Sm,BlockSpr) + f(Sp,BlockSmr));
   BlockSzr = f(Sz,BlockI);
   BlockSpr = f(Sp,BlockI);
   BlockSmr = f(Sm,BlockI);
   
   BlockIold = BlockI;
   BlockI = f(BlockI,m_id2);
   
   HSB = f(BlockSz,BlockSzr)+0.5*(f(BlockSp,BlockSmr)+f(BlockSm,BlockSpr));
   HSB = HSB+f(BlockH,BlockI)+f(BlockI,BlockHr);
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
   
   
   BlockH = T'*BlockH*T;
   BlockSz = T'*BlockSz*T; BlockSp = T'*BlockSp*T; 
   BlockSm = T'*BlockSm*T;
   
   BlockHr = Tr'*BlockHr*Tr;
   BlockSzr = Tr'*BlockSzr*Tr;BlockSpr = Tr'*BlockSpr*Tr; 
   BlockSmr = Tr'*BlockSmr*Tr;
   
   fprintf('%d\t\t %d\t\t\t %f\t\t\t %f\n',i, 2*i+2, Energy/2/(i+1), TruncationError);
end