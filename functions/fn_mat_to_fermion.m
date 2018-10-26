function result = fn_mat_to_fermion(T)
    result=T;
    for alpha1=1:2
        for alpha2=1:2
            for beta1=1:2
                for beta2=1:2
                    if mod((alpha1+beta1)*(alpha2+1),2)==1
                        result(alpha1,alpha2,beta1,beta2) = -1*T(alpha1,alpha2,beta1,beta2);
                    end
                end
            end
        end
    end
% ����������״���Ϊ���ָ�����ż�Ի��һ�����š�  
    
end