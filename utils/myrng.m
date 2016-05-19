function myrng()
% Author		: Tiep Vu (http://www.personal.psu.edu/thv102/)
% Time created	: Wed Jan 27 00:23:58 2016
% Last modified	: Wed Jan 27 00:23:59 2016
% Description	: make sure that random function does not return 
%                 the same output

    c = clock; %get current time 

    t = mod(floor(c(6))*13, 100); % c(6) is 'second'
    for i = 1: t
        randi(t);
    end
end