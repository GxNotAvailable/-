function dX = angCW(t, X, Tp, Te, lambda0)
% ������΢�ֶԲ��й���CW����

% ���볣��
global CONSTANTS
omg = CONSTANTS.refOmega;

% ��ʱ��t�õ���ǰʱ�̵�Э̬��
PhiMtx = lambdaTranMtx(omg, t);
lbd = PhiMtx * lambda0;

% �ɸ�ʱ��Э̬������׷��˫��������
r = norm([lbd(4), lbd(5), lbd(6)]);
if abs(r) < 1e-6
    alpha = 0;
    beta = 0;
else
    beta = asin(lbd(6)/r);
    calpha = lbd(4) / (norm([lbd(4), lbd(5)]));
    salpha = lbd(5) / (norm([lbd(4), lbd(5)]));
    alpha = atan2(salpha, calpha);
end

% CW���̣�������Ϊ���������

% ״̬΢�ַ���
dX = zeros(6,1);
dX(1) = X(4);
dX(2) = X(5);
dX(3) = X(6);
dX(4) = 2*omg*X(5) + 3*omg*omg*X(1) + (Te-Tp) * cos(alpha)*cos(beta);
dX(5) = -2*omg*X(4) + (Te - Tp) * sin(alpha)*cos(beta);
dX(6) = -omg*omg*X(3) + (Te - Tp) * sin(beta);
