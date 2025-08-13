%% ����׷��������������׷�ӹ켣�����ƹ���

%% ��������
traPur = load('CWMinPackPur.txt');
traEva = load('CWMinPackEva.txt');
control = load('CWMinPackControl.txt');
% ׷�������������Ĺ켣
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
% ��Ծ���������ٶ�
DR = tT;
DV = tT;
for i = 1 : length(tT)
    DR(i) = norm([xE(i)-xP(i),yE(i)-yP(i),zE(i)-zP(i)]);
    DV(i) = norm([xVE(i)-xVP(i),yVE(i)-yVP(i),zVE(i)-zVP(i)]);
end
% ׷�������������Ŀ����ɣ������Ƶ�����Ӧ��ͬ
cT = control(:,1);
len = length(cT);
cT = cT(1:len-1,1);
alpha_deg = control(1:len-1,4);
beta_deg = control(1:len-1,5);

%% ����׷����ά�켣
figure;
hold all;
plot3(xP(1)/1000,yP(1)/1000,zP(1)/1000,'sr','LineWidth',3.0);
plot3(xE(1)/1000,yE(1)/1000,zE(1)/1000,'ob','LineWidth',3.0);
legend('׷�ٺ�����','���ݺ�����');
plot3(xP/1000,yP/1000,zP/1000,'-r','LineWidth',3.0);
plot3(xE/1000,yE/1000,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('x/km');
ylabel('y/km');
zlabel('z/km')

%% ���Ƹ�ͶӰƽ����λ�÷����仯ͼ
% xyƽ��
figure;
hold all;
plot(xP(1)/1000,yP(1)/1000,'sr','LineWidth',3.0);
plot(xE(1)/1000,yE(1)/1000,'ob','LineWidth',3.0);
legend('׷�ٺ�����','���ݺ�����');
plot(xP/1000,yP/1000,'-r','LineWidth',3.0);
plot(xE/1000,yE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('x/km');
ylabel('y/km');
% xzƽ��
figure;
hold all;
plot(xP(1)/1000,zP(1)/1000,'sr','LineWidth',3.0);
plot(xE(1)/1000,zE(1)/1000,'ob','LineWidth',3.0);
legend('׷�ٺ�����','���ݺ�����');
plot(xP/1000,zP/1000,'-r','LineWidth',3.0);
plot(xE/1000,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('x/km');
ylabel('z/km');
% yzƽ��
figure;
hold all;
plot(yP(1)/1000,zP(1)/1000,'sr','LineWidth',3.0);
plot(yE(1)/1000,zE(1)/1000,'ob','LineWidth',3.0);
legend('׷�ٺ�����','���ݺ�����');
plot(yP/1000,zP/1000,'-r','LineWidth',3.0);
plot(yE/1000,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('y/km');
ylabel('z/km');

%% ���Ƹ����������ʱ��仯ͼ
% xƽ��
figure;
hold all;
plot(tT,xP/1000,'-r','LineWidth',3.0);
plot(tT,xE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('x/km');
legend('׷�ٺ�����','���ݺ�����','Location','NorthWest');
% yƽ��
figure;
hold all;
plot(tT,yP/1000,'-r','LineWidth',3.0);
plot(tT,yE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('y/km');
legend('׷�ٺ�����','���ݺ�����','Location','SouthEast');
% zƽ��
figure;
hold all;
plot(tT,zP/1000,'-r','LineWidth',3.0);
plot(tT,zE/1000,'-b','LineWidth',3.0);
box on;
grid on;
xlabel('t/s');
ylabel('z/km');
legend('׷�ٺ�����','���ݺ�����');

%% ������Ծ�������ٶȱ仯����
figure;
plot(tT, DR/1000, '-b','LineWidth',3.0);
box on;
ylabel('\DeltaR / km')
figure;
plot(tT, DV, '-b','LineWidth',3.0);
box on;
ylabel('\DeltaV / (m/s)')

%% ���ƿ����ɱ仯
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
