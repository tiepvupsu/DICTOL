function [st] = getTimeStr()
    c = clock();
    MM = num2str(c(2) + 100); MM = MM(2:end);
    dd = num2str(c(3) + 100); dd = dd(2:end);
    hh = num2str(c(4) + 100); hh = hh(2:end);
    mm = num2str(c(5) + 100); mm = mm(2:end);
    ss = num2str(round(c(6)) + 100); ss = ss(2:end);


    st = strcat(MM,dd, '_', hh, mm, ss);
    fprintf('Time now: %s/%s/%4d - %s:%s:%s\n', MM, dd, c(1), hh, mm, ss);
    % if nargin == 2

end 