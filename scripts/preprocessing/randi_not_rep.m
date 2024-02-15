
function idx = randi_not_rep(IMAX,N)
idx = zeros(N,1);
indices = 1:IMAX;
% s = rng;
% rng(seed);
for i = 1:N
    ii = (randi(length(indices),1,1));
    idx(i) = indices(ii);
    indices(ii) = [];
end
% rng(s);