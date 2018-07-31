function [area, perimeter] = portCalculator( port)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

switch port.type
    case "circ" %cylinder configuration
        %diameter refers to the combustion chamber diameter
        area=pi*(port.diameter(end)/2)^2;
        perimeter = pi*port.diameter(end);
               
    case "Dport" %dport configuration

        %arguments contains D_outer, fuel web t, and metal thickness tau
        [area,perimeter] = DPort(port.D_outer,port.fuelweb(end),port.tau);
end

end