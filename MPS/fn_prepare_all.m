function m_mps=fn_prepare_all(mps)
% ����һ��cell��ÿһ�����ǹ̶�������������mps
N=length(mps);

m_mps=cell(1,N);

for k=1:N
    
    m_mps{k}=fn_prepare_othersitemps(mps,k);

end

end