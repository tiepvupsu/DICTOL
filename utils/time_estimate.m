function time_estimate(t)
    t = round(t);
    h = floor(t/3600);
    t = t - 3600*h;
    m = floor(t/60);
    t = t - m*60;
    fprintf('| time left: %2dh%2dm%2ds\n', h, m, t);
end 