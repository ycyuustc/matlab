
v_lambda = v_lambda_store;
v_y = v_result_store;

v_alpha0 = [3,2.6,1.5,8];
fn_error(v_alpha0,v_lambda,v_y);

v_alpha = v_alpha0;
error_now = fn_error(v_alpha,v_lambda,v_y);
step = 0.1;
tv = eye(4);
while 1
     
     v_flag = zeros(1,4);
    for i=1:4
       
        while 1
            v_delta_alpha = tv(i,:)*step;
            v_alpha_try_left = v_alpha - v_delta_alpha;
            v_alpha_try_right = v_alpha + v_delta_alpha;
            error_try_left = fn_error(v_alpha_try_left,v_lambda,v_y);
            error_try_right = fn_error(v_alpha_try_right,v_lambda,v_y);
            if error_try_left<error_now
                v_alpha = v_alpha_try_left;
                error_now = error_try_left;
                disp(error_now);
                disp(v_alpha);
                v_flag(i) = v_flag(i)+1;
            else
                if error_try_right<error_now
                    v_alpha = v_alpha_try_right;
                    error_now = error_try_right;
                    disp(error_now);
                    disp(v_alpha);
                     v_flag(i) = v_flag(i)+1;
                else
                    break;
                end
            end
            
        end
    end
    
    if sum(v_flag) == 0
        step = step/2;
    end
    
    disp(['step=', num2str(step)]);
    if step<1e-6
        break;
    end
    
end

disp(v_alpha);
[~,v_r] = fn_error(v_alpha,v_lambda,v_y);
close(figure(2));
figure(2);
hold on;
plot(v_lambda,v_y,'-r*');
plot(v_lambda,v_r,'-bo');





function [res,v_r] = fn_error(v_alpha,v_lambda,v_y)

alpha1 = v_alpha(1)*1e5;
alpha2 = v_alpha(2);
alpha3 = v_alpha(3)*1e6;
alpha4 = v_alpha(4);

Num = 2001;
vx = linspace(650,730,Num);
v_length = vx(2:end)-vx(1:end-1);
vx_nt = zeros(1,Num);
vx_nt(1:Num-1) = 1/2*v_length;
vx_nt(2:Num) = vx_nt(2:Num) + 1/2*v_length;

Num_lambda = length(v_lambda);
m_x = ones(Num_lambda,1)*vx;
m_lambda = v_lambda'*ones(1,Num);

s=1.51;
% alpha1 = 3e5;
% alpha2 = 2.6;
% alpha3 = 1.5e6;
% alpha4 = 8;

m_n = alpha1./m_lambda.^2 + alpha2;
m_alpha = power(10,alpha3./m_lambda.^2-alpha4);

m_A = 16.*m_n.^2.*s;
m_B = (m_n+1).^3.*(m_n+s^2);
m_C = 2.*(m_n.^2-1).*(m_n.^2-s^2);
m_D = (m_n-1).^3.*(m_n-s^2);

m_ex = exp(-m_alpha.*m_x);
m_integral = (m_A.*m_ex)./...
    (m_B-m_C.*cos(4*pi.*m_n.*m_x./m_lambda).*m_ex+m_D.*m_ex.^2);

v_r = sum(m_integral.*vx_nt,2)/80;

res = sum((v_r-v_y').^2);

end

