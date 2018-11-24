function m_result=fn_contract_perm1(mps)

N=length(mps);

m_result=zeros(N,N);
norm=fn_contractmps(mps);

for m=1:N-1
for n=m+1:N
   
    %%%%%%%%%% Part I: generate L  %%%%%%%%%%    
    % step 1 :
    L=fn_contract(conj(mps{1}),3,1,mps{1},3,1);
    L=permute(L,[1,3,2,4]);
    % step 2 :
    for k=2:m
       L=fn_contract(L,4,1,conj(mps{k}),3,1);
       L=fn_contract(L,5,[1,5],mps{k},3,[1,3]);
       L=permute(L,[3,4,1,2]);
    end
    
    %%%%%%%%%% Part II: generate L+ %%%%%%%%%%
    Lplus=fn_contract(L,4,[1,4],conj(mps{m+1}),3,[1,3]);
    Lplus=permute(Lplus,[3,1,2]);
    
    %%%%%%%%%% Part III: generate R %%%%%%%%%%
    % step 1: 
    R=1;
    % step 2:
    for k=N:-1:(n+1)
       R=fn_contract(conj(mps{k}),3,2,R,2,1);
       R=fn_contract(mps{k},3,[2,3],R,3,[3,2]);
       R=permute(R,[2,1]);
    end
    
    %%%%%%%%%% Part IV: genreate R+ %%%%%%%%%%
    Rplus=fn_contract(mps{n},3,2,R,2,2);
    Rplus=permute(Rplus,[3,1,2]);
    
    %%%%%%%%%% Part V: expand L+ %%%%%%%%%%
    for k=(m+1):(n-1)
        Lplus=fn_contract(Lplus,3,1,conj(mps{k+1}),3,1);
        Lplus=fn_contract(Lplus,4,[1,4],mps{k},3,[1,3]);
        Lplus=permute(Lplus,[2,3,1]);
    end
    
    %%%%%%%%%% Part IV: last step, contract L+ and R+ %%%%%%%%%%
    result=fn_contract(Lplus,3,[1,2,3],Rplus,3,[1,2,3]);
    
    m_result(m,n)=result;
    m_result(n,m)=conj(result);
end
end

m_result(logical(eye(N)))=norm;
   
end
