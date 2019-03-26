lieb1 = [2,1,4,3,0];

lieb2 = [2,1,4,3,0];
lieb0 = fn_lieb_prod(lieb1,lieb2);
disp(lieb0);
disp(closed_circle);

disp(fn_lieb_contract(lieb0));
% disp(fn_lieb_add(lieb1,3,2));
M = 6;
lieb0 = [(M+1):(2*M),1:M,0];
tm = [5,6;1,2;3,4;5,6;3,4;2,3;4,5];
[num,~] = size(tm);
for i = 1:num
    p1 = tm(i,1);
    p2 = tm(i,2);
    lieb0 = fn_lieb_add(lieb0,p1,p2);
end
disp(lieb0);


N = 4;
lambda = 1;
v_ind = fn_generate_poisson(lambda,N);
lieb = [(N+1):(2*N),1:N,0];
for i=1:length(v_ind)
    lieb = fn_lieb_add_base(lieb,v_ind(i));
end
disp(lieb);


beta = 6;
M = 6;
N = 8;
lambda = beta/M;

c_lieb = cell(1,M);
storage_up = cell(1,M);
storage_down = cell(1,M);
lieb_identity = fn_lieb_identity(N);

% initialization
for i=1:M
    c_lieb{i} = lieb_identity;
    storage_up{i} = lieb_identity;
    storage_down{i} = lieb_identity;
end
circle_number = N;
lieb_tot = lieb_identity;
% end initialization

fn_virsualization(c_lieb,1);
test_lieb =cell(1,M);
for m = 1:M
    test_lieb{m} = fn_generate_lieb(1,N);
end

fn_virsualization(test_lieb,2);
tot_prod_lieb = test_lieb{1};
for m = 2:M
    tot_prod_lieb = fn_lieb_prod(tot_prod_lieb,test_lieb{m});
end
fn_virsualization({tot_prod_lieb},3);




% begin the random walk
Number_step = 10000;
m_store = zeros(Number_step,2*N+1);
index_step = 1;
while 1
    
    for index = 1:(M-1)
        
        lieb_random = fn_generate_lieb(lambda,N);
        lieb_up = storage_up{index};
        lieb_down = storage_down{index};
        lieb_tot_try = fn_lieb_prod(lieb_up,lieb_random);
        lieb_tot_try = fn_lieb_prod(lieb_tot_try,lieb_down);
        circle_number_try = fn_lieb_contract(lieb_tot_try);
        prob = 1/(1+2^(circle_number-circle_number_try));
        
        if rand()<prob
            circle_number = circle_number_try;
            lieb_tot = lieb_tot_try;
            c_lieb{index} = lieb_random;
        end
        storage_up{index+1} = ...
            fn_lieb_prod(storage_up{index},c_lieb{index});
        
        m_store(index_step,:) = lieb_tot;
        index_step = index_step +1;
        if index_step>Number_step
            break;
        end
    end
    if index_step>Number_step
        break;
    end
    
    for index = M:-1:2
        
        lieb_random = fn_generate_lieb(lambda,N);
        lieb_up = storage_up{index};
        lieb_down = storage_down{index};
        lieb_tot_try = fn_lieb_prod(lieb_up,lieb_random);
        lieb_tot_try = fn_lieb_prod(lieb_tot_try,lieb_down);
        circle_number_try = fn_lieb_contract(lieb_tot_try);
        prob = 1/(1+2^(circle_number-circle_number_try));
        
        if rand()<prob
            circle_number = circle_number_try;
            lieb_tot = lieb_tot_try;
            c_lieb{index} = lieb_random;
        end
        storage_down{index-1} = ...
            fn_lieb_prod(c_lieb{index},storage_down{index});
        
        m_store(index_step,:) = lieb_tot;
        index_step = index_step +1;
        if index_step>Number_step
            break;
        end
    end
    if index_step>Number_step
        break;
    end
    
    disp(index_step);
    
end  % walking circle end

disp('haha');


%%%%%%%%   END  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         the functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function res = fn_lieb_identity(N)
res = [(N+1):(2*N),1:N,0];
end

function res_lieb = fn_lieb_prod(lieb1,lieb2)

closed_circle = lieb1(end)+lieb2(end);
M = (length(lieb1)-1)/2;
lieb1 = lieb1(1:(2*M));
lieb2 = lieb2(1:(2*M));
lieb0 = zeros(1,2*M);

v_ind = find(lieb1(1:M)<=M);
lieb0(v_ind) = lieb1(v_ind);
v_ind = find(lieb2((M+1):(2*M))>=M+1);
lieb0(M+v_ind) = lieb2(M+v_ind);

v_judge = zeros(1,M);
point = 1;
while point<=M
    if lieb0(point)==0
        p = point;
        while 1
            p = lieb1(p);
            if p<=M
                break;
            end
            p = p - M;
            v_judge(p) = 1;
            p = lieb2(p);
            if p>=M+1
                break;
            end
            v_judge(p) = 1;
            p = p + M;
        end
        lieb0(point) = p;
        lieb0(p) = point;
    end
    
    point = point + 1;
end

point = M+1;
while point<=2*M
    if lieb0(point)==0
        p = point;
        while 1
            p = lieb2(p);
            if p>M
                break;
            end
            v_judge(p) = 1;
            p = p + M;
            p = lieb1(p);
            if p<=M
                break;
            end
            p = p - M;
            v_judge(p) = 1;
        end
        lieb0(point) = p;
        lieb0(p) = point;
    end
    
    point = point + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% counting the circle number
point = 1;
while point<=M
    if v_judge(point)==0
        p = point;
        while 1
            p = p + M;
            p = lieb1(p);
            p = p - M;
            v_judge(p) = 1;
            p = lieb2(p);
            v_judge(p) = 1;
            if p==point
                break;
            end
        end
        closed_circle = closed_circle + 1;
    end
    
    point = point + 1;
end
res_lieb = [lieb0,closed_circle];
end

function res_lieb=fn_lieb_add(lieb0,p1,p2)

M = (length(lieb0)-1)/2;
closed_circle = lieb0(end);
lieb0 = lieb0(1:(2*M));

q1 = lieb0(p1 + M);
q2 = lieb0(p2 + M);

if q2==p1+M
    lieb = lieb0;
    closed_circle = closed_circle + 1;
else
    lieb = lieb0;
    lieb(q1) = q2;
    lieb(q2) = q1;
    lieb(p1+M) = p2+M;
    lieb(p2+M) = p1+M;
    %     closed_circle = closed_circle + 0;
end
res_lieb = [lieb,closed_circle];
end

function lieb = fn_lieb_add_base(lieb0,position)
M = (length(lieb0)-1)/2;
p1 = position;
p2 = mod(position,M) + 1;
lieb = fn_lieb_add(lieb0,p1,p2);
end

function circle_number = fn_lieb_contract(lieb0)

circle_number = lieb0(end);
M = (length(lieb0)-1)/2;
lieb0 = lieb0(1:(2*M));
v_finished = zeros(1,2*M);

point = 1;
while point<=2*M
    
    if v_finished(point)==0
        p = point;
        while 1
            p = lieb0(p);
            v_finished(p) = 1;
            if p<=M
                p = p + M;
            else
                p = p - M;
            end
            v_finished(p) = 1;
            if p==point
                break;
            end
        end
        circle_number = circle_number + 1;
    end
    
    point = point + 1;
end

end

function res_ind = fn_generate_poisson(lambda,N)
%poisson distribution
num_leg = random('Poisson',lambda*N);
if num_leg>0
    v_rand = rand(1,num_leg)*N;
    v_ind = floor(v_rand+1);
    v_position = mod(v_rand,1);
    [~,v_sort] = sort(v_position);
    v_ind = v_ind(v_sort);
    res_ind = v_ind;
else
    res_ind = [];
end

end

function res_lieb = fn_generate_lieb(lambda,N)

v_ind = fn_generate_poisson(lambda,N);
% disp(v_ind);
if isempty(v_ind)
    res_lieb = fn_lieb_identity(N);
else
    lieb = fn_lieb_identity(N);
    for i=1:length(v_ind)
        lieb = fn_lieb_add_base(lieb,v_ind(i));
        res_lieb = lieb;
    end
end

end

function fn_virsualization(c_lieb,figure_number)
close(figure(figure_number));
figure(figure_number);
hold on;

M = length(c_lieb);
N = (length(c_lieb{1})-1)/2;

figure(figure_number);
axis([0,N+1,0,M+2]);

for y = 1:(M+1)
    vx = linspace(0.5,(N+0.5),2);
    vy = y*ones(1,2);
    plot(vx,vy,'r--');
    for x=1:N
        figure(figure_number);
        plot(x,y,'ro');
    end
end

for n=1:N
    x1 = n; y1 = 1;
    x2 = n; y2 = 1/2;
    plot([x1,x2],[y1,y2],'-b','linewidth',3);
    x1 = n; y1 = M+1;
    x2 = n; y2 = M+3/2;
    plot([x1,x2],[y1,y2],'-b','linewidth',3);
end

for m=1:M
    y_center = M + 3/2 - m;
    lieb = c_lieb{m};
    lieb = lieb(1:end-1);
    v_judge = zeros(1,2*N);
    point = 1;
    while point<=2*N
        if v_judge(point)==0
            p1 = point;
            p2 = lieb(p1);
            v_judge(p1) = 1;
            v_judge(p2) = 1;
            
            if p1<=N
                x1 = p1;
                y1 = y_center + 1/2;
            else
                x1 = p1 - N;
                y1 = y_center - 1/2;
            end
            
            if p2<=N
                x2 = p2;
                y2 = y_center + 1/2;
            else
                x2 = p2 - N;
                y2 = y_center - 1/2;
            end
            
            if (p1>N && p2<=N) || (p1<=N && p2>N)
                plot([x1,x2],[y1,y2],'-b','linewidth',3);
            elseif p1<=N
                max_x = max(x1,x2);
                min_x = min(x1,x2);
                A = (max_x-min_x)/N*0.5;
                num_x = 100;
                vx = linspace(min_x,max_x,num_x);
                vy = (y1+y2)/2 ...
                    - A*sin((vx-min_x)./(max_x-min_x)*pi);
                plot(vx,vy,'-b','linewidth',3);
            else
                max_x = max(x1,x2);
                min_x = min(x1,x2);
                A = (max_x-min_x)/N*0.5;
                num_x = 100;
                vx = linspace(min_x,max_x,num_x);
                vy = (y1+y2)/2 ...
                    + A*sin((vx-min_x)./(max_x-min_x)*pi);
                plot(vx,vy,'-b','linewidth',3);
            end
        end
        
        point = point + 1;
    end
end

end

