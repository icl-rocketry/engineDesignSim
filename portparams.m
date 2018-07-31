
function [area, perimeter] = portparams(porttype, arguments)

switch porttype
    case 1 %cylinder configuration
        D_outer=arguments(1);
        fuelweb=arguments(2);
        
        a=D_outer/2-fuelweb; %radius of combustion chamber port
        
        area = pi*a^2;
        perimeter = 2*pi*a;
        
        
    case 2 %dport configuration
        %arguments contains D_outer, fuel web t, and metal thickness tau
        [area,perimeter] = DPort(arguments(1),arguments(2),arguments(3));
end

end
