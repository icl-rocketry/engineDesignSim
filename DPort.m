function [area, perimeter] = DPort(D_outer, t, tau)

%D_outer = 60;
%t = 20;
a = D_outer/2-t;
%tau = 2;

theta = acos((t+tau)/a);

area = 2*(a^2 * theta - a *(t+tau)*sin(theta));

perimeter = 4*a*(theta + sin(theta));

end


