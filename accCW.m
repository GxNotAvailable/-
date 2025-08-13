function dX = accCW(t, X, tVec, axVec, ayVec, azVec)
% �̶�������΢�ֶԲ��й���CW����

% ���볣��
global CONSTANTS
omg = CONSTANTS.refOmega;

% CW���̣�������ֱ��Ϊ�������ٶȣ����ò�ֵ���
ux =  interp1(tVec,axVec,t,'spline'); 
uy =  interp1(tVec,ayVec,t,'spline'); 
uz =  interp1(tVec,azVec,t,'spline'); 

% ״̬΢�ַ���
dX = zeros(6,1);
dX(1) = X(4);
dX(2) = X(5);
dX(3) = X(6);
dX(4) = 2*omg*X(5) + 3*omg*omg*X(1) + ux;
dX(5) = -2*omg*X(4) + uy;
dX(6) = -omg*omg*X(3) + uz;
