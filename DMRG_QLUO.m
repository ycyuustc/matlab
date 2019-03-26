%% %%%%%%%%%%%%%%%%%%%%%%%%%%% WARN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This MATLAB code is written By Qiang Luo@RUC.
% It is mainly for personal use only.
% For more information, contact:
%   Q. Luo  - qiangluo@ruc.edu.cn
% DATA & PALACE: Apr. 23 2016, at CSRC.
%% %%%%%%%%%%%%%%%%%%%% INPUT PARAMETER(S) %%%%%%%%%%%%%%%%%%%%%%%
% clc,clear
NIter = 10;             % 迭代次数，系统的链长则 Ltot = 2*NIter+2
m = 20;                 % 保留态数目，测试情况下一般不超过50.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Intialize local operators 
Sz = [1/2 0;0 -1/2];    %sigma_z
Sp = [0 1;0 0];         %s^+， Sp=Sx+Im*Sy;   
Sm = [0 0;1 0];         %s^-， Sm=Sx-Im*Sy;  
dim=2;
I= eye(dim);
Hmid = kron(Sz,Sz)+0.5*(kron(Sp,Sm)+kron(Sm,Sp));
% 构建初始时系统块和环境块的哈密顿量和自旋算符
BlockH  = zeros(dim);
BlockSz = Sz;BlockSp = Sp;BlockSm = Sm;

BlockHr  = zeros(dim);
BlockSzr = Sz;BlockSpr = Sp;BlockSmr = Sm;

NKeep = dim;
fprintf('Iter\tTotalLength\t Energy\t\t\t\t Trunc\n');
%====================【Begin main iterations】=============================
for l = 1:NIter
    BlockI = eye(NKeep);
    %%  step 1:  构建迭代时系统块和环境块的哈密顿量和自旋算符
    BlockH = kron(BlockH, I) + kron(BlockSz, Sz) + 0.5 * ( kron(BlockSp, Sm) + kron(BlockSm, Sp) );
    BlockSz = kron(BlockI, Sz);BlockSp = kron(BlockI, Sp);BlockSm = kron(BlockI, Sm);
    
    BlockHr = kron(I, BlockHr) + kron(Sz, BlockSzr) + 0.5 * ( kron(Sm, BlockSpr) + kron(Sp, BlockSmr) );
    BlockSzr = kron(Sz, BlockI);BlockSpr = kron(Sp, BlockI);BlockSmr = kron(Sm, BlockI);
    
    BlockIold = BlockI;
    BlockI  = kron(BlockI, I);
    %%  step 2:  HAMILTONIAN MATRIX for superblock
    HSB = kron(BlockSz, BlockSzr) + 0.5 * ( kron(BlockSp, BlockSmr) + kron(BlockSm, BlockSpr) );
    HSB = HSB + kron(BlockH, BlockI) + kron(BlockI, BlockHr);
    HSB = 0.5 * (HSB + HSB');  % ensure H is symmetric 
    opts.disp  = 0; opts.issym = 1; opts.real  = 1;
    [Psi, Energy] = eigs(HSB,1,'SA', opts);  % Diagonalizing the Hamiltonian, which is a real, symmetric matrix.
  
    %  Form the reduced density matrix
    len=length(Psi);
    PsiMatrix = reshape(Psi, sqrt(len), sqrt(len));
    PsiMatrix = transpose(PsiMatrix);        % matlab是列主元的，装置之后行（i）和列（j）分别对应于系统块和环境块的基矢。
    
    Rho = PsiMatrix * PsiMatrix';
    [V, D] = eig(Rho);
    [D, Index] = sort(diag(D), 'descend');   % sort eigenvalues descending
    V = V(:,Index);                          % sort eigenvectors the same way
    NKeep = min(size(D, 1), m);
    T = V(:, 1:NKeep);

    Rhor = PsiMatrix' * PsiMatrix;
    [Vr, Dr] = eig(Rhor);
    [Dr, Indexr] = sort(diag(Dr), 'descend');  % sort eigenvalues descending
    Vr = Vr(:,Indexr);                         % sort eigenvectors the same way
    NKeepr = min(size(Dr, 1), m);
    Tr = Vr(:, 1:NKeepr);
    
    TruncationError = 1 - sum(D(1:NKeep));    
    %%  step 3: Transform the block operators into the truncated basis
    BlockH  = T'*BlockH*T;
    BlockSz = T'*BlockSz*T; BlockSp = T'*BlockSp*T; BlockSm = T'*BlockSm*T;
    
    BlockHr  = Tr'*BlockHr*Tr;
    BlockSzr = Tr'*BlockSzr*Tr; BlockSpr = Tr'*BlockSpr*Tr; BlockSmr = Tr'*BlockSmr*Tr;
    
    fprintf('%d\t\t %d\t\t\t %f\t\t\t %f\n',l, 2*l+2, Energy/(2*l+2), TruncationError);
end