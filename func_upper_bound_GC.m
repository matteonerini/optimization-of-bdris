function P = func_upper_bound_GC(hIT, hRI, NG)
% The received power upper bound is returned

NI = size(hIT,1);

if NG == 1 % Single connected
    P = sum(abs(hRI' .* hIT)) ^ 2;   
elseif NG == 0 % Fully connected
    P = norm(hRI) ^ 2 * norm(hIT) ^ 2;
else % Group connected
    try
        hRI_norm = vecnorm(reshape(hRI,[NG,NI/NG]));
        hIT_norm = vecnorm(reshape(hIT,[NG,NI/NG]))';
        P = (hRI_norm * hIT_norm) ^ 2;
    catch
        P = nan;
    end
end
end