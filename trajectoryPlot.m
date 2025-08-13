%% 绘制追踪器与逃逸器的追逃轨迹、控制轨线

%% 读入数据
traPur = load('CWMinPackPur.txt');
traEva = load('CWMinPackEva.txt');
control = load('CWMinPackControl.txt');
% 追踪器与逃逸器的轨迹
tT = traPur(:,1);
xP = traPur(:,2);
yP = traPur(:,3);
zP = traPur(:,4);
xVP = traPur(:,5);
yVP = traPur(:,6);
zVP = traPur(:,7);
xE = traEva(:,2);
yE = traEva(:,3);
zE = traEva(:,4);
xVE = traEva(:,5);
yVE = traEva(:,6);
zVE = traEva(:,7);
% 相对距离与相对速度
DR = tT;
DV = tT;
for i = 1 : length(tT)
    DR(i) = norm([xE(i)-xP(i),yE(i)-yP(i),zE(i)-zP(i)]);
    DV(i) = norm([xVE(i)-xVP(i),yVE(i)-yVP(i),zVE(i)-zVP(i)]);
end
% 追踪器与逃逸器的控制律，根据推导它们应相同
cT = control(:,1);
len = length(cT);
cT = cT(1:len-1,1);
alpha_deg = control(1:len-1,4);
beta_deg = control(1:len-1,5);

%% 绘制追逃三维轨迹
figure;
hold all;
plot3(xP(1)/1000,yP(1)/1000,zP(1)/1000,'sr','LineWidth',3.0);
plot3(xE(1)/1000,yE(1)/1000,zE(1)/1000,'ob','LineWidth',3.0);
legend('追踪航天器','逃逸航天器');
plot3(xP/1000,yP/1000,zP/1000,'-r','LineWidth',3.0);
plot3(xE/1000,yE/1000,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('x/km');
ylabel('y/km');
zlabel('z/km')

%% 绘制各投影平面内位置分量变化图
% xy平面
figure;
hold all;
plot(xP(1)/1000,yP(1)/1000,'sr','LineWidth',3.0);
plot(xE(1)/1000,yE(1)/1000,'ob','LineWidth',3.0);
legend('追踪航天器','逃逸航天器');
plot(xP/1000,yP/1000,'-r','LineWidth',3.0);
plot(xE/1000,yE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('x/km');
ylabel('y/km');
% xz平面
figure;
hold all;
plot(xP(1)/1000,zP(1)/1000,'sr','LineWidth',3.0);
plot(xE(1)/1000,zE(1)/1000,'ob','LineWidth',3.0);
legend('追踪航天器','逃逸航天器');
plot(xP/1000,zP/1000,'-r','LineWidth',3.0);
plot(xE/1000,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('x/km');
ylabel('z/km');
% yz平面
figure;
hold all;
plot(yP(1)/1000,zP(1)/1000,'sr','LineWidth',3.0);
plot(yE(1)/1000,zE(1)/1000,'ob','LineWidth',3.0);
legend('追踪航天器','逃逸航天器');
plot(yP/1000,zP/1000,'-r','LineWidth',3.0);
plot(yE/1000,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('y/km');
ylabel('z/km');

%% 绘制各坐标分量随时间变化图
% x平面
figure;
hold all;
plot(tT,xP/1000,'-r','LineWidth',3.0);
plot(tT,xE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('x/km');
legend('追踪航天器','逃逸航天器','Location','NorthWest');
% y平面
figure;
hold all;
plot(tT,yP/1000,'-r','LineWidth',3.0);
plot(tT,yE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('y/km');
legend('追踪航天器','逃逸航天器','Location','SouthEast');
% z平面
figure;
hold all;
plot(tT,zP/1000,'-r','LineWidth',3.0);
plot(tT,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('z/km');
legend('追踪航天器','逃逸航天器');

%% 绘制相对距离相对速度变化规律
figure;
plot(tT, DR/1000, '-b','LineWidth',3.0);
box on;
ylabel('\DeltaR / km')
figure;
plot(tT, DV, '-b','LineWidth',3.0);
box on;
ylabel('\DeltaV / (m/s)')

%% 绘制控制律变化
figure;
plot(cT,alpha_deg,'-r','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('\alpha/deg');
figure;
plot(cT,beta_deg,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('\beta/deg');
