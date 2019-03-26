disp(fn_g([1,1]));

Num = 501;
vx1 = linspace(-0.5,0.5,Num);
vx2 = linspace(-0.5,0.5,Num);
[x1,x2] = meshgrid(vx1,vx2);

delta_x = 1/(Num-1);
m_nt = ones(Num,Num);
m_nt(1,1) = 0.25;
m_nt(1,Num) = 0.25;
m_nt(Num,1) = 0.25;
m_nt(Num,Num) = 0.25;
m_nt(1,2:(Num-1)) = 0.5;
m_nt(Num,2:(Num-1)) = 0.5;
m_nt(2:(Num-1),1) = 0.5;
m_nt(2:(Num-1),Num) = 0.5;

m_nt = m_nt * delta_x^2;
k1 = rand()*4;
k2 = k1;
f = exp(1j*(k1*x1+k2*x2)).*(sign(x2-x1)+1)/2;

res_int_1 = sum(sum(f.*m_nt));
res_int_2 = exp(-1j/2*(k1+k2))*fn_g([k1+k2,k2,0]);
disp(res_int_1);
disp(res_int_2);

k3 = rand();
disp(fn_g([k1+k2+k3,k2+k3,k3,0]));
disp(fn_g([k1+k2+k3,k2+k3,k2,0]));


m_permutation = perms([3,2,1]);
v_sigma = [1,-1];
vx = [1,2,3];
vk1 = kk_t12(:,1);
vk2 = kk_t12(:,2);

v_momentum = linspace(-25,25,1000);
v_res = zeros(1,1000);

for im = 1:1000
    momentum = v_momentum(im);  
    my_sum = 0;
    for ip = 1:6       
        v_sort = m_permutation(ip,:);        
        for i1 = 1:3
            for i2 = 1:3
                for si11 = 1:2
                    for si21 = 1:2
                        for si12 = 1:2
                            for si22 = 1:2
                                
                                vx(v_sort) = [1,2,3];
                                if vx(3)>vx(2)
                                    r1 = 1;
                                else
                                    r1 = 2;
                                end
                                
                                if vx(1)>vx(2)
                                    r2 = 1;
                                else
                                    r2 = 2;
                                end
                                
                                sigma11 = v_sigma(si11);
                                sigma12 = v_sigma(si12);
                                sigma21 = v_sigma(si21);
                                sigma22 = v_sigma(si22);
                                
                                k1 = sigma12*vk1(i2) - momentum;
                                k2 = -sigma21*vk2(i1) + sigma22*vk2(i2);
                                k3 = -sigma11*vk1(i1) + momentum;
                                vk = [k1,k2,k3];
                                vk = vk(v_sort);
                                k1 = vk(1); k2 = vk(2); k3 = vk(3);
                                
                                my_sum = my_sum + ...
                                    fn_g([k1+k2+k3,k2+k3,k3,0])*...
                                    A(r2,i2,si12,si22)*...
                                    conj(A(r1,i1,si11,si21))*...
                                    exp(-1j/2*(k1+k2+k3));
                                
                            end
                        end
                    end
                end
            end
        end
        
    end
    
    v_res(im) = my_sum;
    
    disp(im);disp(my_sum);
    
end

figure(1);
plot(v_momentum,real(v_res),'bo');




function res = fn_g(vk)

if length(vk)>1
    [kmax,max_ind] = max(vk);
    [kmin,min_ind] = min(vk);
    
    if kmax-kmin<1e-6
        res = exp(1j*mean(vk));
    else
        vk1 = vk; vk2 = vk;
        vk1(min_ind) = []; vk2(max_ind) = [];
        res = (fn_g(vk1)-fn_g(vk2))/1j/(kmax-kmin);
    end
    
else
    res = exp(1j*vk(1));
end

end