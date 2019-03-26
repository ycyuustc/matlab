% m_conf = ones(8,8);

M = 10;
N = 8;

m_conf = diag(power(-1,0:(M-1)))*ones(M,N);

m_conf(7,2) = -m_conf(7,2);
% m_conf(6,2) = -m_conf(6,2);
% m_conf(4,3) = -m_conf(4,3);
% m_conf(2,4) = -m_conf(2,4);

M = 4; N = 4;

m_conf = diag(power(-1,0:(M-1)))*ones(M,N);
fn_visualization(m_conf,2);
m_conf(4,2) = -m_conf(4,2);
% m_conf(3,3) = -m_conf(3,3);
% m_conf(2,4) = -m_conf(2,4);
% m_conf(6,2) = -m_conf(6,2);
% m_conf(4,3) = -m_conf(4,3);
% m_conf(2,4) = -m_conf(2,4);

fn_visualization(m_conf,3);

lieb = [5,6,7,8,1,2,3,4,0];
fn_lieb_visualization({lieb},4);

lieb0 = fn_lieb_identity(6);
lieb1 = fn_lieb_add_base(lieb0,4);
lieb2 = fn_lieb_add_base(lieb0,2);
lieb3 = fn_lieb_add_base(lieb0,2);
lieb4 = fn_lieb_add_base(lieb0,1);
lieb5 = fn_lieb_add_base(lieb0,5);

fn_lieb_visualization({lieb1,lieb2,lieb3,lieb4,lieb5},5);

function fn_visualization(m_conf,plotnum)
close(figure(plotnum));
figure(plotnum);
hold on;

[M,N] = size(m_conf);
figure(plotnum); axis([0,N+1,0,M+1]);
for m = 1:M
    for n = 1:N
        
        figure(plotnum);
        
        x0 = n;
        y0 = M+1-m;
        d1 = 1/4; d2 = 1/2;
        x1 = x0 - d1; y1 = y0; v1 = [x1,y1];
        x2 = x0; y2 = y0 + d1; v2 = [x2,y2];
        x3 = x0 + d1; y3 = y0; v3 = [x3,y3];
        x4 = x0; y4 = y0 - d1; v4 = [x4,y4];
        x5 = x0 - d2; y5 = y0; v5 = [x5,y5];
        x6 = x0; y6 = y0 + d2; v6 = [x6,y6];
        x7 = x0 + d2; y7 = y0; v7 = [x7,y7];
        x8 = x0; y8 = y0 - d2; v8 = [x8,y8];
        
        if m_conf(m,n) == -1
            fn_point_plot(v5,v1,v2,v6);
            fn_point_plot(v8,v4,v3,v7);
        end
        
        if m_conf(m,n) == 1
           fn_point_plot(v5,v1,v4,v8);
           fn_point_plot(v6,v2,v3,v7);
        end
        
        plot(x0,y0,'ro');
        
    end
end

end



function fn_point_plot(v1,v2,v3,v4)

x1 = v1(1); y1 = v1(2);
x2 = v2(1); y2 = v2(2);
x3 = v3(1); y3 = v3(2);
x4 = v4(1); y4 = v4(2);

plot([x1,x2],[y1,y2],'b','linewidth',3);
plot([x2,x3],[y2,y3],'b','linewidth',3);
plot([x3,x4],[y3,y4],'b','linewidth',3);

end



function fn_lieb_visualization(c_lieb,figure_number)
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

function [res_lieb,edg_number] = fn_generate_lieb(lambda,N)

v_ind = fn_generate_poisson(lambda,N);
edg_number = length(v_ind);
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
