function dP = mRiccati(t, P, A, S, Q)
% �����迨��΢�ַ���
P = reshape(P, size(A));
dP = -A'*P - P*A + P*S*P - Q;
dP = dP(:);
