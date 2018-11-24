function [v_1b,m_2b]=fn_datamining_mps(mps)

N=length(mps);
v_1b=cell(1,N);
m_2b=cell(N,N);

for k=1:N
    tensor=mps{k};
    L=fn_contractnetL(mps,k);
    R=fn_contractnetR(mps,k);
    L=fn_contract(L,2,1,conj(tensor),3,1);
    L=fn_contract(L,3,1,tensor,3,1);
    LR=fn_contract(L,4,[1,3],R,2,[1,2]);
    v_1b{k}=LR;
    m_2b{k,k}=LR;
end

for i=1:N-1
    for j=i+1:N
   
        L=fn_contractnetL(mps,i);
        R=fn_contractnetR(mps,j);
        L=fn_contract(L,2,1,conj(mps{i}),3,1);
        L=fn_contract(L,3,1,mps{i},3,1);
        R=fn_contract(conj(mps{j}),3,2,R,2,1);
        R=fn_contract(R,3,3,mps{j},3,2);
        
        k=i+1;
        while k<j-0.1
           L=fn_contract(L,4,1,conj(mps{k}),3,1);
           L=fn_contract(L,5,[2,5],mps{k},3,[1,3]);
           L=permute(L,[3,1,4,2]);
           k=k+1;
        end
        
        tensor=fn_contract(L,4,[1,3],R,4,[1,3]);
        m_2b{i,j}=tensor;
        m_2b{j,i}=permute(tensor,[3,4,1,2]);
        
    end
end

