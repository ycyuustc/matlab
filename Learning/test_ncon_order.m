%input variables
beta = 0.6;
max_level = 6;
chi_u = 6;
chi_w = 6;
chi_v = 6;
chi_y = 6;
UL_max = 100;
loop_max = 100;

%convergence criteria
UL_convergence = 1e-10;
B_convergence = 1e-7;
D_convergence = 1e-7;
A_convergence = 1e-7;

%build Ising transfer matrix. Break the symmetry with small perturbation
epsilon = 0.1;
Q = [(1+epsilon/4)*exp(beta) exp(-beta); exp(-beta) (1-epsilon/4)*exp(beta)];

%define delta
delta3 = zeros(2,2,2);
for i=1:2
    delta3(i,i,i) = 1;
end

%create plaquette operator
A = ncon({delta3,Q,delta3,Q,delta3,Q,delta3,Q},...
    {[-1,8,1],[1,2],[-2,2,3],[3,4],[-3,4,5],[5,6],[-4,6,7],[7,8]});

%normalise
A = A./sqrt(ncon({A,conj(A)},{[1,2,3,4],[1,2,3,4]}));

u = random('unif',-1,1,[min(chi_u,size(A,1)),size(A,1),min(chi_u,size(A,1)),size(A,1)]);
    w = random('unif',-1,1,[min(chi_w,size(u,3)*size(u,1)),size(u,3),size(u,1)]);
    v_l = random('unif',-1,1,[min(chi_v,size(A,4)*size(u,1)),size(A,4),size(u,1)]);
    v_r = random('unif',-1,1,[min(chi_v,size(A,2)*size(u,1)),size(u,3),size(A,2)]);
    y_l = random('unif',-1,1,[min(chi_y,size(v_l,1)*size(v_l,1)),size(v_l,1),size(v_l,1)]);
    y_r = random('unif',-1,1,[min(chi_y,size(v_r,1)*size(v_r,1)),size(v_r,1),size(v_r,1)]);


B_old = mcon({v_l,u,v_r,A,A,conj(A),conj(A),conj(u),conj(v_l),conj(v_r)},...
    {[-1,1,2],[2,5,3,6],[-2,3,4],[5,9,7,1],[6,4,8,9],[11,12,7,10],[13,14,8,12],[15,11,16,13],[-4,10,15],[-3,16,14]});

function res = mcon(cTensor,cIndex)

numTen = length(cIndex);
numIndex = 0;
numIndexOut = 0;
for i=1:numTen

    numIndex = numIndex + length(cIndex{i});
    numIndexOut = numIndexOut + sum(cIndex{i}<0); 
    disp(numIndexOut);
    disp(numIndex);
    
end

numIndexIn = (numIndex - numIndexOut)/2;

v_order(netcon(cIndex,0,2,1,1)) = 1:numIndexIn;
c_order = cell(1,numTen);
for i=1:numTen
    
   tv = cIndex{i};
   for j = 1:length(tv)
       if tv(j)>0
           tv(j) = v_order(tv(j));
       end
   end
   
   c_order{i} = tv;
   
end

res = ncon(cTensor,c_order);

end