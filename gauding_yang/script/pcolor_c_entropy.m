Num = 600;

v_c = linspace(0.8, 1.25, Num);
v_T0 = linspace(0.5e-5,4e-3,Num);

[m_c,m_T0] = meshgrid(v_c,v_T0);

m_r = zeros(Num,Num);
h0 = 0.5;
mu0 = -0.2;
for k1=1:Num
    for k2=1:Num
        
        coupling = m_c(k1,k2);
        T0 = m_T0(k1,k2);
        
        h = h0/coupling^2;
        mu = mu0/coupling^2;
        T = T0/coupling^2;
        
        m_r(k1,k2) = pl_gaudin_entropy(h,mu,T)*coupling;
        
    end
    disp(k1);
    
end