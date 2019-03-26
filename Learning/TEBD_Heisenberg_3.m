% ABS->ASB algorithm

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

epsilon = -0.01;
mh = -1/2*(0.0*hzz + 1.0*hz);
h12 = reshape(mh,[2,2,2,2]);
hm = expm(epsilon*mh);
hamil = reshape(hm,[2,2,2,2]);

D =2;
% initialization
A = rand(2,2,2);
B = rand(2,2,2);
S = diag(rand(1,2));
% end initialization

tot_H = 10;
num_step = floor(tot_H/abs(epsilon));
for kkk = 1:num_step
    % update A-B
    [A,S,B] = fn_ABS(A,B,S,hamil,D);
    % finished updating A-B
    
    % update B-A
    [B,S,A] = fn_ABS(B,A,S,hamil,D);
    % finished updating B-A
    
%     disp(diag(sAB));
%     disp(diag(sBA));
    
    vAB = diag(S).^2;
    entropy_AB_pre = entropy_AB;
    entropy_AB = -sum(vAB.*log(vAB));
    disp(['entropy_AB = ',num2str(entropy_AB)]);
    
    % calculate the energy
    TAB = ncon({A,B,S},{[-1,1,-3],[1,2,-4],[2,-2]});
    TBA = ncon({B,S,A},{[-1,1,-3],[1,2],[2,-2,-4]});
    
    T2AB = ncon({TAB,conj(TAB)},{[-1,-3,1,2],[-2,-4,1,2]});
    T2BA = ncon({TBA,conj(TBA)},{[-1,-3,1,2],[-2,-4,1,2]});
    ts2 = size(T2AB); ts = [ts2,ones(1,4-length(ts2))];
    mAB = reshape(T2AB,ts(1)*ts(2),ts(3)*ts(4));
    [v,~] = eigs(mAB,1);
    mAB_inf = v*v';
    TAB_inf = reshape(mAB_inf,ts);
    TAB_local = mcon({TAB,conj(TAB),TAB_inf},...
        {[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
    ts2 = size(T2BA); ts = [ts2,ones(1,4-length(ts2))];
    mBA = reshape(T2BA,ts(1)*ts(2),ts(3)*ts(4));
    [v,~] = eigs(mBA,1);
    mBA_inf = v*v';
    TBA_inf = reshape(mBA_inf,ts);
    TBA_local = mcon({TBA,conj(TBA),TBA_inf},...
        {[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
       
    energy_AB = ncon({TAB_local,h12},{[1,2,3,4],[1,2,3,4]})/...
        ncon({TAB_local},{[1,2,1,2]});
    disp(energy_AB);
    energy_BA = ncon({TBA_local,h12},{[1,2,3,4],[1,2,3,4]})/...
        ncon({TBA_local},{[1,2,1,2]});
    disp(energy_BA);
    
    plot(kkk,energy_AB,'ro');
    plot(kkk,energy_BA,'bo');
    
    pause(0.001);
    %
end


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

function [A,S,B] = fn_ABS(A,B,S,U,D)

T = ncon({A,B,S,U},{[-1,1,3],[1,2,4],[2,-3],[-2,-4,3,4]});
ts = size(T);
mT = reshape(T,ts(1)*ts(2),ts(3)*ts(4));
[u,s,v] = svd(mT,'econ'); v = v';
truncate_D = min([D,sum(diag(s)/norm(diag(s))>1e-6)]);
mA = u(:,1:truncate_D);
mS = s(1:truncate_D,1:truncate_D);
mB = v(1:truncate_D,:);
A = reshape(mA,[ts(1),ts(2),truncate_D]);
A = permute(A,[1,3,2]);
S = mS/norm(diag(mS));
B = reshape(mB,[truncate_D,ts(3),ts(4)]);

end

