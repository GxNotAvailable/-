% ������΢�ֶԲ���⣬�ؼ��ǹ��������ֵ���Ⲣ�������
clear all
clc

%% ��������
% ������
earthRadius         =   6378.137 * 1000;                % ����뾶��m 
gravParam           =   3.986004418E14;                 % ��������������m3/s2
% �ο�Բ�������
refHeight           =   36000 * 1000;                     % �ο�����߶ȣ�m
refRadius           =   earthRadius + refHeight;        % �ο�����볤�ᣬm
refOmega            =   sqrt(gravParam / (refRadius*refRadius*refRadius)); % �ο�������ٶȣ�rad/s
% ׷�ٺ����������ݺ�������ʼ�������
Xp = [-6*1000; -16*1000; 4*1000; -9; 13.6; 0];
Xe = [0; 0; 0; 0; 0; 0];
% ׷������������������������ٶ�
g0 = 9.78;                      % ��ƽ���������ٶ�
TpMax = 0.06;
TeMax = 0.04;
% �Բ߳�ʼ״̬����
X0 = Xe - Xp;
% ΢�ֶԲ߻��ֲ���
stepT               =   10;                              % s
% ����������ʼ����ٶȺͳ�ʼ��Ծ���
r0 = norm([X0(1), X0(2), X0(3)]);
v0 = norm([X0(4), X0(5), X0(6)]);
% �򵥹��Ƶõ�׷��ʱ��
a0 = TpMax - TeMax;
estT = (sqrt(v0*v0+2*a0*r0) - v0) / a0;

%% ����λ
% ��������λ����Ҫ����
scales.length       =   1;%refRadius;                      % ��������λ���ο�����뾶
scales.speed        =   1;%sqrt(gravParam/scales.length);  % �ٶ�����λ����������λ��Ӧ��Բ����ٶȴ�С
% ��������λ���ɻ�������λ����
scales.time         =   scales.length/scales.speed;     % ʱ������λ������/�ٶ�
scales.acceleration =   scales.speed/scales.time;       % ���ٶ�����λ���ٶ�/ʱ��
scales.omega        =   1/scales.time;                  % ���ٶ�����λ��1/ʱ��
scales.gravparam    =   scales.acceleration*scales.length^2;    % ������������λ�����ٶ�*����^2
% ����֮��ĳ�����ȫ�ֱ������������ļ���ʹ��
global CONSTANTS
CONSTANTS.refOmega  =   refOmega / scales.omega;        % ����֮��Ĳο�������ٶ�
CONSTANTS.TpMax = TpMax / scales.acceleration;
CONSTANTS.TeMax = TeMax / scales.acceleration;          % ����֮���������������������ٶ�
CONSTANTS.stepT     =   stepT / scales.time;            % ����֮��Ļ��ֲ���
CONSTANTS.X0 = X0;
CONSTANTS.X0(1:3) = X0(1:3) / scales.length;
CONSTANTS.X0(4:6) = X0(4:6) / scales.speed;
CONSTANTS.Xp = Xp;
CONSTANTS.Xp(1:3) = Xp(1:3) / scales.length;
CONSTANTS.Xp(4:6) = Xp(4:6) / scales.speed;
CONSTANTS.Xe = Xe;
CONSTANTS.Xe(1:3) = Xe(1:3) / scales.length;
CONSTANTS.Xe(4:6) = Xe(4:6) / scales.speed;
CONSTANTS.estT = estT / scales.time;

%% ������߱�ֵ���⣬�����������һ�������Է�����
%% ���ý����㷨��ó�ֵ
% ����GA


% ����DE


% ����PSO
npars = 4;          % ��Ʊ�����4��
lb = [600/scales.time -10 -10 -10];
ub = [1200/scales.time 10 10 10];
optPSO.npop = npars*15;         optPSO.niter = 300;
optPSO.Display = 1;          optPSO.Plot = 0;
optPSO.Vectorized = 0;     optPSO.Parallel = 0;
optPSO.StallIterLimit = 100;
[yopt, fval, exitflag, output] = ppso(@psoFunc,npars,lb,ub, optPSO);

% % �Ա���Ϊ׷��ʱ��tf���Լ�Э̬����ֵ��v1,v2,v3����Ҫ������ֵ
Y0 = yopt;
options = optimoptions('fsolve','Display','iter', 'TolFun', 1e-10, 'TolX', 1e-10, 'MaxFunEvals', 1000); % Option to display output
[Y,fval] = fsolve(@endGame, Y0, options) % Call solver




