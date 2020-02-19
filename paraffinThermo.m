function [T_flame, gamma, m_mol, R,c_star] = paraffinThermo(etac)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

T_flame = 2740;

gamma = 1.2;

m_mol = 0.028;

R = 8.314/m_mol;

%c_star = 1406;

c_star = etac*sqrt(gamma*R*T_flame)/(gamma*(2/(gamma+1))^((gamma+1)/(2*gamma-2))); %spad

end

