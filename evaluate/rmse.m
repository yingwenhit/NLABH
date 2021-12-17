function s = rmse(u, Io)
[M, N] = size(u);
MN = M*N;
s = sqrt((sum(sum((u-Io).^2)))/(MN));
end