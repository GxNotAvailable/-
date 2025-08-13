function Mtx = lambdaTranMtx(omg, deltaT)
% 构造协态变量的状态转移矩阵

theta = omg * deltaT;
s = sin(omg * deltaT);
c = cos(omg * deltaT);

Mtx = zeros(6,6);

Mtx(1,1) = 4 - 3*c;
Mtx(1,2) = 6*theta - 6*s;
Mtx(1,4) = -3*omg*s;
Mtx(1,5) = -6*omg + 6*omg*c;

Mtx(2,2) = 1;

Mtx(3,3) = c;
Mtx(3,6) = omg*s;

Mtx(4,1) = -s/omg;
Mtx(4,2) = 2*(c-1)/omg;
Mtx(4,4) = c;
Mtx(4,5) = 2*s;

Mtx(5,1) = 2*(1-c)/omg;
Mtx(5,2) = (3*theta-4*s)/omg;
Mtx(5,4) = -2*s;
Mtx(5,5) = 4*c - 3;

Mtx(6,3) = -s/omg;
Mtx(6,6) = c;
