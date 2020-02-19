function [area, perimeter] = portCalculator( port)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

area=pi*(port.diameter(end)/2)^2;
perimeter = pi*port.diameter(end);

end