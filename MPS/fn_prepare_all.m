function m_mps=fn_prepare_all(mps)
% 返回一个cell，每一个都是固定点的正交化后的mps
N=length(mps);

m_mps=cell(1,N);

for k=1:N
    
    m_mps{k}=fn_prepare_othersitemps(mps,k);

end

end