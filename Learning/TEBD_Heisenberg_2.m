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
mh = -1/2*(hzz + 0.00*hx + 1.01/2*hz);
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
entropy_AB = 0;
entropy_BA = 0;

tot_H = 30;
num_step = floor(tot_H/abs(epsilon));
for kkk = 1:num_step
    % update A-B
    TATB = ncon({TA,sAB,TB},{[-1,1,-2],[1,2],[2,-3,-4]});
    TATB = ncon({TATB,hamil},{[-1,1,-3,2],[-2,-4,1,2]});
    ts = size(TATB);
    matrix = reshape(TATB,ts(1)*ts(2),ts(3)*ts(4));
%     disp(max(max(abs(matrix))));
    [u,s,v] = svd(matrix,'econ'); v = v';
    truncate_D = min(D,length(diag(s)));
    mA = u(:,1:truncate_D);
    mS = s(1:truncate_D,1:truncate_D);
    mB = v(1:truncate_D,:);
    TA = reshape(mA,[ts(1),ts(2),truncate_D]);
    TA = permute(TA,[1,3,2]);
    TB = reshape(mB,[truncate_D,ts(3),ts(4)]);
    sAB = mS;
    sAB = sAB/norm(diag(sAB));
    
    [sBA,TA,sAB,TB] = fn_TEBD_canonical_2_v2(sBA,TA,sAB,TB);
    tt = ncon({sBA,TA,sAB,TB,sBA},{[-1,1],[1,2,-3],[2,3],[3,4,-4],[4,-2]});
    rho = ncon({tt,conj(tt)},{[1,2,-3,-4],[1,2,-1,-2]});
    disp('energy1');
    disp(ncon({rho,h12},{[1,2,3,4],[1,2,3,4]})/ncon({rho},{[1,2,1,2]}));
    % finished updating A-B
    
    % update B-A
    TBTA = ncon({TB,sBA,TA},{[-1,1,-2],[1,2],[2,-3,-4]});
    TBTA = ncon({TBTA,hamil},{[-1,1,-3,2],[-2,-4,1,2]});
    ts = size(TBTA);
    matrix = reshape(TBTA,ts(1)*ts(2),ts(3)*ts(4));
%     disp(max(max(abs(matrix))));
    [u,s,v] = svd(matrix,'econ'); v = v';
    truncate_D = min(D,length(diag(s)));
    mA = u(:,1:truncate_D);
    mS = s(1:truncate_D,1:truncate_D);
    mB = v(1:truncate_D,:);
    TB = reshape(mA,[ts(1),ts(2),truncate_D]);
    TB = permute(TB,[1,3,2]);
    TA = reshape(mB,[truncate_D,ts(3),ts(4)]);
    sBA = mS;
    sBA = sBA/norm(diag(sBA));
    
    [sAB,TB,sBA,TA] = fn_TEBD_canonical_2_v2(sAB,TB,sBA,TA);
    tt = ncon({sAB,TB,sBA,TA,sAB},{[-1,1],[1,2,-3],[2,3],[3,4,-4],[4,-2]});
    rho = ncon({tt,conj(tt)},{[1,2,-3,-4],[1,2,-1,-2]});
    disp('energy2');
    disp(ncon({rho,h12},{[1,2,3,4],[1,2,3,4]})/ncon({rho},{[1,2,1,2]}));
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
    ts = size(T2);
    m2 = reshape(T2,ts(1)*ts(2),ts(3)*ts(4));
    [v,d] = eigs(m2,1);
    m_inf = v*v';
    T_inf = reshape(m_inf,ts);
    T_local = mcon({T,conj(T),T_inf},{[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
    energy = ncon({T_local,h12},{[1,2,3,4],[1,2,3,4]})/...
        ncon({T_local},{[1,2,1,2]});
    disp('energy energy 1');
    disp(energy);
    plot(kkk/num_step*tot_H,energy,'ro');
%     correlation0 = ncon({T_local,thzz},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
%     correlation1 = ncon({T_local,hz1},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
%     correlation2 = ncon({T_local,hz2},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
%     disp(['correlation = ',num2str(correlation0-correlation1*correlation2)]);
%     
    T = ncon({TB,sBA,TA,sAB},{[-1,1,-3],[1,2],[2,3,-4],[3,-2]});
    T2 = ncon({T,conj(T)},{[-1,-3,1,2],[-2,-4,1,2]});
    ts = size(T2);
    m2 = reshape(T2,ts(1)*ts(2),ts(3)*ts(4));
    [v,d] = eigs(m2,1);
    m_inf = v*v';
    T_inf = reshape(m_inf,ts);
    T_local = mcon({T,conj(T),T_inf},{[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
    energy = ncon({T_local,h12},{[1,2,3,4],[1,2,3,4]})/...
        ncon({T_local},{[1,2,1,2]});
    disp('energy energy 2'); 
    disp(energy);
    plot(kkk/num_step*tot_H,energy,'bo');
    pause(0.001);
%     correlation0 = ncon({T_local,thzz},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
%     correlation1 = ncon({T_local,hz1},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
%     correlation2 = ncon({T_local,hz2},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
%     disp(['correlation = ',num2str(correlation0-correlation1*correlation2)]);
    %
end


T1 = ncon({TA,sAB},{[-1,1,-3],[1,-2]});
T2 = ncon({TB,sBA},{[-1,1,-3],[1,-2]});



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

