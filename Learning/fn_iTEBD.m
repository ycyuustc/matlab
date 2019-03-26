function [TA,sAB,TB,sBA] = fn_iTEBD(TA,sAB,TB,sBA,h12,D,delta,num_walk)

close(figure(258));
figure(258);
hold on;

ts = size(h12);
mh12 = reshape(h12,ts(1)*ts(2),ts(3)*ts(4));
mhamil = expm(-delta*mh12);
hamil = reshape(mhamil,ts);


for kkk = 1:num_walk
    
   % update A-B
    TATB = ncon({TA,sAB,TB},{[-1,1,-2],[1,2],[2,-3,-4]});
    TATB = ncon({TATB,hamil},{[-1,1,-3,2],[-2,-4,1,2]});
    ts = size(TATB);
    matrix = reshape(TATB,ts(1)*ts(2),ts(3)*ts(4));
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
    energy1=ncon({rho,h12},{[1,2,3,4],[1,2,3,4]})/ncon({rho},{[1,2,1,2]});
    disp(['energy1=',num2str(energy1)]);
    % finished updating A-B
    
    % update B-A
    TBTA = ncon({TB,sBA,TA},{[-1,1,-2],[1,2],[2,-3,-4]});
    TBTA = ncon({TBTA,hamil},{[-1,1,-3,2],[-2,-4,1,2]});
    ts = size(TBTA);
    matrix = reshape(TBTA,ts(1)*ts(2),ts(3)*ts(4));
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
    energy2=ncon({rho,h12},{[1,2,3,4],[1,2,3,4]})/ncon({rho},{[1,2,1,2]});
    disp(['energy2=',num2str(energy2)]);
    % finished updating B-A
    
    vAB = diag(sAB).^2;
    vBA = diag(sBA).^2;
    entropy_AB = -sum(vAB.*log(vAB));
    entropy_BA = -sum(vBA.*log(vBA));
    disp(['entropy_AB = ',num2str(entropy_AB)]);
    disp(['entropy_BA = ',num2str(entropy_BA)]);
    
    % calculate the energy
%     T = ncon({TA,sAB,TB,sBA},{[-1,1,-3],[1,2],[2,3,-4],[3,-2]});
%     T2 = ncon({T,conj(T)},{[-1,-3,1,2],[-2,-4,1,2]});
%     ts = size(T2);
%     m2 = reshape(T2,ts(1)*ts(2),ts(3)*ts(4));
%     [v,~] = eigs(m2,1);
%     m_inf = v*v';
%     T_inf = reshape(m_inf,ts);
%     T_local = mcon({T,conj(T),T_inf},{[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
%     energy = ncon({T_local,h12},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
%     disp(energy);
%     plot(kkk,energy,'ro');
% 
%     T = ncon({TB,sBA,TA,sAB},{[-1,1,-3],[1,2],[2,3,-4],[3,-2]});
%     T2 = ncon({T,conj(T)},{[-1,-3,1,2],[-2,-4,1,2]});
%     ts = size(T2);
%     m2 = reshape(T2,ts(1)*ts(2),ts(3)*ts(4));
%     [v,~] = eigs(m2,1);
%     m_inf = v*v';
%     T_inf = reshape(m_inf,ts);
%     T_local = mcon({T,conj(T),T_inf},{[2,4,-3,-4],[1,3,-1,-2],[4,3,2,1]});
%     energy = ncon({T_local,h12},{[1,2,3,4],[1,2,3,4]})/...
%         ncon({T_local},{[1,2,1,2]});
     
%     disp(energy);
    plot(kkk,energy1,'ro');
    plot(kkk,energy2,'b+');
    pause(0.0001);
    
end

end