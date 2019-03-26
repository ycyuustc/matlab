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
% mh = -1/2*(hxx + 1.01/2*hz);
mh = hzz + 0.01*hz;
h12 = reshape(mh,[2,2,2,2]);
hm = expm(epsilon*mh);
hamil = reshape(hm,[2,2,2,2]);

D =64;
% initialization
TA = rand(2,2,2);
TB = rand(2,2,2);
sAB = diag(rand(1,2));
sBA = diag(rand(1,2));
% end initialization
entropy_AB = 0;
entropy_BA = 0;

tot_H = 30;
num_step = floor(tot_H/abs(epsilon));
for kkk = 1:num_step
    % update A-B
    TATB = ncon({sBA,TA,sAB,TB,sBA},{[-1,1],[1,2,-3],[2,3],[3,4,-4],[4,-2]});
    temp = ncon({TATB,hamil},{[-1,-2,1,2],[-3,-4,1,2]});
    TATB = ncon({TATB,hamil},{[-1,-3,1,2],[-2,-4,1,2]});
    ts = size(TATB);
    matrix = reshape(TATB,ts(1)*ts(2),ts(3)*ts(4));
%     disp(max(max(abs(matrix))));
    [u,s,v] = svd(matrix,'econ'); v = v';
    truncate_D = min([D,sum(diag(s)/norm(diag(s))>1e-6)]);
    mA = u(:,1:truncate_D);
    mS = s(1:truncate_D,1:truncate_D);
    mB = v(1:truncate_D,:);
    TA = reshape(mA,[ts(1),ts(2),truncate_D]);
    TA = permute(TA,[1,3,2]);
    TB = reshape(mB,[truncate_D,ts(3),ts(4)]);
    sBA_inv = diag(1./diag(sBA));
    TA = ncon({sBA_inv,TA},{[-1,1],[1,-2,-3]});
    TB = ncon({TB,sBA_inv},{[-1,1,-3],[1,-2]});
    sAB = mS;
    temp_a = norm(diag(sAB));
    sAB = sAB/norm(diag(sAB));
    % finished updating A-B
    temp2 = ncon({sBA,TA,sAB,TB,sBA},{[-1,1],[1,2,-3],[2,3],[3,4,-4],[4,-2]});
    temp2 = temp2*temp_a;
    % update B-A
    TBTA = ncon({sAB,TB,sBA,TA,sAB},{[-1,1],[1,2,-3],[2,3],[3,4,-4],[4,-2]});
    temp = ncon({TBTA,hamil},{[-1,-2,1,2],[-3,-4,1,2]});
    TBTA = ncon({TBTA,hamil},{[-1,-3,1,2],[-2,-4,1,2]});
    ts = size(TBTA);
    matrix = reshape(TBTA,ts(1)*ts(2),ts(3)*ts(4));
%     disp(max(max(abs(matrix))));
    [u,s,v] = svd(matrix,'econ'); v = v';
    truncate_D = min([D,sum(diag(s)/norm(s)>1e-6)]);
    mA = u(:,1:truncate_D);
    mS = s(1:truncate_D,1:truncate_D);
    mB = v(1:truncate_D,:);
    TB = reshape(mA,[ts(1),ts(2),truncate_D]);
    TB = permute(TB,[1,3,2]);
    TA = reshape(mB,[truncate_D,ts(3),ts(4)]);
    sAB_inv = diag(1./diag(sAB));
    TB = ncon({sAB_inv,TB},{[-1,1],[1,-2,-3]});
    TA = ncon({TA,sAB_inv},{[-1,1,-3],[1,-2]});
    sBA = mS;
    temp_a = norm(diag(sBA));
    sBA = sBA/norm(diag(sBA));
    temp2 = ncon({sAB,TB,sBA,TA,sAB},{[-1,1],[1,2,-3],[2,3],[3,4,-4],[4,-2]});
    temp2 = temp2*temp_a;
    % finished updating B-A
    
%     disp(diag(sAB));
%     disp(diag(sBA));
    
    vAB = diag(sAB).^2;
    vBA = diag(sBA).^2;
    entropy_AB_pre = entropy_AB;
    entropy_BA_pre = entropy_BA;
    entropy_AB = -sum(vAB.*log(vAB));
    entropy_BA = -sum(vBA.*log(vBA));
    disp(['entropy_AB = ',num2str(entropy_AB)]);
    disp(['entropy_BA = ',num2str(entropy_BA)]);
    
%     figure(1);hold on;
%     plot(vAB,'ro');
%     if abs(entropy_AB_pre-entropy_AB)<1e-8 ...
%             && abs(entropy_BA_pre-entropy_BA)<1e-8
%         break;
%     end   

    % calculate the energy
    T = ncon({TA,sAB,TB,sBA},{[-1,1,-3],[1,2],[2,3,-4],[3,-2]});
    T2 = ncon({T,conj(T)},{[-1,-3,1,2],[-2,-4,1,2]});
    ts2 = size(T2); ts = [ts2,ones(1,4-length(ts2))];
    m2 = reshape(T2,ts(1)*ts(2),ts(3)*ts(4));
    [v,~] = eigs(m2,1);
    m_inf = v*v';
    T_inf = reshape(m_inf,ts);
    T_local = mcon({T,conj(T),T_inf},{[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
    energy = ncon({T_local,h12},{[1,2,3,4],[1,2,3,4]})/...
        ncon({T_local},{[1,2,1,2]});
    disp(energy);
    plot(kkk,energy,'ro');
    
    T = ncon({TB,sBA,TA,sAB},{[-1,1,-3],[1,2],[2,3,-4],[3,-2]});
    T2 = ncon({T,conj(T)},{[-1,-3,1,2],[-2,-4,1,2]});
    ts2 = size(T2); ts = [ts2,ones(1,4-length(ts2))];
    m2 = reshape(T2,ts(1)*ts(2),ts(3)*ts(4));
    [v,~] = eigs(m2,1);
    m_inf = v*v';
    T_inf = reshape(m_inf,ts);
    T_local = mcon({T,conj(T),T_inf},{[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
    energy = ncon({T_local,h12},{[1,2,3,4],[1,2,3,4]})/...
        ncon({T_local},{[1,2,1,2]});
     
    disp(energy);
    plot(kkk,energy,'bo');
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

