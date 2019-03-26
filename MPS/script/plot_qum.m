tm1 = -m_r(:,:,1);
tm2 = -m_r(:,:,2);
tm3 = -m_r(:,:,3);
tm4 = -m_r(:,:,4);
tm5 = -m_r(:,:,5);
tm6 = -m_r(:,:,6);

close(figure(1));
figure(1);
hold on;
plot(vN,tm1(8,:),'--*');
plot(vN,tm2(8,:),'--*');
plot(vN,tm3(8,:),'--*');
plot(vN,tm4(8,:),'--*');
plot(vN,tm5(8,:),'--*');
plot(vN,tm6(8,:),'--*');

close(figure(2));
figure(2);
hold on;
plot(vN,tm1(8,:),'-*');
plot(v_N,m_eig(:,1),'ro');

close(figure(3));
figure(3);
hold on;
plot(vN,tm2(8,:),'-*');
plot(v_N,m_eig(:,2),'ro');

close(figure(4));
figure(4);
hold on;
plot(vN,tm3(8,:),'-*');
plot(v_N,m_eig(:,3),'ro');

close(figure(5));
figure(5);
hold on;
plot(vN,tm4(8,:),'-*');
plot(v_N,m_eig(:,4),'ro');

close(figure(6));
figure(6);
hold on;
plot(vN,tm5(8,:),'-*');
plot(v_N,m_eig(:,5),'ro');

close(figure(7));
figure(7);
hold on;
plot(vN,tm6(8,:),'-*');
plot(v_N,m_eig(:,6),'ro');
