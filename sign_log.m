function log_sign_m=sign_log(m)
    log_sign_m = zeros(size(m));
    for i = 1:size(m, 1)
        for j = 1:size(m, 3)
            if m(i, j) > 0
                log_sign_m(i, j) = log(m(i, j));
            elseif m(i, j) < 0
                log_sign_m(i, j) = log(-m(i, j));
            end
        end
    end
end