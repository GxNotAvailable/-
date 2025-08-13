function F = endGame(y)
% ������΢�ֶԲ��й������߱�ֵ����������
% �Ա���Ϊ׷��ʱ��tf��Э̬����ֵ�е�v1,v2,v3
% ��ʽԼ��Ϊ�ն�λ��Լ�����ն˹��ܶ����������Լ��

global CONSTANTS
omg = CONSTANTS.refOmega;

% �Ա����е�׷�ӽ���ʱ���Э̬����ֵ
tf = y(1);
lbd_f = zeros(6,1);
lbd_f(1) = y(2);
lbd_f(2) = y(3);
lbd_f(3) = y(4);

% ��Э̬����ֵ���Э̬����ʼֵ
PhiMtx_t02tf = lambdaTranMtx(omg, tf);      % ��ʼ���ն�Э̬��״̬ת�ƾ���
if abs(det(PhiMtx_t02tf)) < 1e-10
    F = zeros(4,1);
    F(1) = 1e4;
    F(2) = 1e4;
    F(3) = 1e4;
    F(4) = 1e4;
    return;
else
    lbd_0 = inv(PhiMtx_t02tf) * lbd_f;               % t0ʱ�̵�Э̬��
    % ���������õ���Э̬����ʼֵ�����ֻ�öԲ�״̬��ֵ
    tspan = [0: CONSTANTS.stepT: tf];
    options_ode = odeset('RelTol', 1e-6, 'AbsTol', 1e-6);
    [TSPAN, XSPAN] = ode45(@(t,X)angCW(t, X, CONSTANTS.TpMax, CONSTANTS.TeMax, lbd_0), tspan, CONSTANTS.X0, options_ode);
    Xf = XSPAN(end,:);

    % �����ն˹��ܶ���
    Hf = lbd_f(1)*Xf(4) + lbd_f(2)*Xf(5) + lbd_f(3)*Xf(6);

    % �����ն˵�ʽԼ��
    F = zeros(4,1);
    F(1) = Xf(1);
    F(2) = Xf(2);
    F(3) = Xf(3);
    F(4) = Hf + 1;
end