function dP = mRiccati(t, P, A, S, Q)
% 矩阵黎卡提微分方程
P = reshape(P, size(A));
dP = -A'*P - P*A + P*S*P - Q;
dP = dP(:);
