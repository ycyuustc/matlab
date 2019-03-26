function [res_cell,res_matrix] = fn_generate_QTM_mps(tau,lambda,Num)

b = @(lambda) 1j/(lambda+1j);
c = @(lambda) lambda/(lambda+1j);

% m2id=eye(2);
% m2z=[1,0;0,-1];
% dT = 2;

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

lambda_1 = tau + lambda;
lambda_2 = tau - lambda;
T1 = c(lambda_1)*m4p + b(lambda_1)*m4id;
T2 = c(lambda_2)*m4p2 + b(lambda_2)*m4id;

res_cell = cell(1,Num);
for i=1:Num
   
    if mod(i,2)==1
       res_cell{i} = T1; 
    else
       res_cell{i} = T2;
    end
    
end

end
