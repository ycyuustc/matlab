% do not contract the sBA at the first step

close(figure(1));
figure(1);
hold on;

h12 = fn_generate_local_heisenberg(2,2,1);
N = 10;
hm = fn_generate_local_heisenberg(N,1,2);
for k = 2:(N-1)
   hm = hm + fn_generate_local_heisenberg(N,k,k+1); 
end
hm = hm + fn_generate_local_heisenberg(N,N,1);
matrix_h = reshape(hm,2^N,2^N);
vh = eig(matrix_h);
energy = min(vh);
energy = energy/N;
disp(energy);

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

%
m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;
mp = reshape(m4p,4,4);
epsilon = -0.001;
hm = expm(epsilon*mp);
hamil = reshape(hm,[2,2,2,2]);
%

sx = [0,1;1,0];
sI = eye(2);
sz = [1,0;0,-1];
hxx = kron(sx,sx);
thxx = reshape(hxx,[2,2,2,2]);
hzz = kron(sz,sz);
hz = kron(sz,sI) + kron(sI,sz);
hx = kron(sx,sI) + kron(sI,sx);
hz1 = kron(sz,sI);
hz2 = kron(sI,sz);
thz1 = reshape(hz1,[2,2,2,2]);
thz2 = reshape(hz2,[2,2,2,2]);
thzz = reshape(hzz,[2,2,2,2]);

epsilon = -0.025;
mh = -1/2*(hxx + 0.00*hx + 1.01/2*hz);
h12 = reshape(mh,[2,2,2,2]);
hm = expm(epsilon*mh);
hamil = reshape(hm,[2,2,2,2]);

D = 8;
% initialization
TA = rand(2,2,2);
TB = rand(2,2,2);
sAB = diag(rand(1,2));
sBA = diag(rand(1,2));
% end initialization

mps = cell(1,4);
mps{1} = TA;
mps{2} = sAB;
mps{3} = TB;
mps{4} = sBA;

mps = fn_iTEBD(mps,h12,8,0.01,10000);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function res = fn_generate_local_heisenberg(N,ind_i,ind_j)

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;
tm = eye(2^N);
c_iden = reshape(tm,2*(ones(1,2*N)));
v_ind = -(1:(2*N));
v_ind(ind_i+N) = 1;
v_ind(ind_j+N) = 2;
v_ind2 = [1,2,-ind_i-N,-ind_j-N];
res = ncon({c_iden,m4p},{v_ind,v_ind2});

end

