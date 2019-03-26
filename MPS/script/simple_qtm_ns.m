m2id=eye(2);
m2z=[1,0;0,-1];

m4id=zeros(2,2,2,2);
m4id(1,1,1,1)=1;
m4id(1,2,1,2)=1;
m4id(2,1,2,1)=1;
m4id(2,2,2,2)=1;

m4p=zeros(2,2,2,2);
m4p(1,1,1,1)=1;
m4p(1,2,2,1)=1;
m4p(2,1,1,2)=1;
m4p(2,2,2,2)=1;

m4p2 = permute(m4p,[1,4,3,2]);

% D=18;
precision = 1e-8;
% N = 10;

% vD = 2:2:16;
% vN = 2:2:20;
% vbeta = [0.1,0.2,0.3,0.4,0.5,0.6];
vD = 20;
vN = 10;
vbeta = 0.1;

m_r = zeros(length(vD),length(vN),length(vbeta));
% m_r2 = zeros(length(vD),length(vN));
for ibeta = 1:length(vbeta)
    beta = vbeta(ibeta);    
    for iN = 1:length(vN)
        for iD = 1:length(vD)
            
            N = vN(iN);
            D = vD(iD);
            
            lambda = beta/N;
            
            p1 = 1/(1-lambda);
            p2 = lambda/(lambda-1);
            
            T1 = p1*m4p + p2*m4id;
            T2 = p1*m4p2 + p2*m4id;
            
            hset=cell(2,N);
            
            T1_left = T1(1,:,:,:);
            T1_left = permute(T1_left,[1,3,2,4]);
            T2_right = T2(:,:,1,:);
            T2_right = permute(T2_right,[1,3,2,4]);
            hset{1,1} = T1_left;
            hset{1,N} = T2_right;
            
            T1_left = T1(2,:,:,:);
            T1_left = permute(T1_left,[1,3,2,4]);
            T2_right = T2(:,:,2,:);
            T2_right = permute(T2_right,[1,3,2,4]);
            hset{2,1} = T1_left;
            hset{2,N} = T2_right;
            
            for i=2:(N-1)
                
                if mod(i,2) == 1
                    hset{1,i} = 1*permute(T1,[1,3,2,4]);
                    hset{2,i} = 1*permute(T1,[1,3,2,4]);
                else
                    hset{1,i} = permute(T2,[1,3,2,4]);
                    hset{2,i} = permute(T2,[1,3,2,4]);
                end
                
            end
            
            hset{1,1} = hset{1,1};
            hset{2,1} = hset{2,1};
            
            [E0,mps0]=fn_minimizeE(hset,D,precision);
            [E0,mps0]=fn_minimizeE_ns(hset,D,precision);
            
            m_r(iD,iN,ibeta) = real(E0);
            
            % [E1,mps1]=fn_minimizeEP(hset,D,precision,mps0);
        end
    end
end

plot_qum;