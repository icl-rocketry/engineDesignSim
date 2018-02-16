function [T_flame, gamma, m_mol, R,c_star] = thermochem(OF,P_cc,etac)

%interp uses X is P_cc and Y is OF

load propepinterp

%T_flame = 3300;
%gamma = 1.24;
%m_mol = 0.0262109;  %Molar mass (kg/mol) should be determined properly


T_flame = interp2(P_cc_vals,OF_vals,T_flame_data,P_cc,OF); %[K] 
gamma = interp2(P_cc_vals,OF_vals,gamma_data,P_cc,OF); %dimensionless
m_mol = interp2(P_cc_vals,OF_vals,m_mol_data,P_cc,OF); %g/mol
%convert m_mol to kg/mol
m_mol=m_mol/1000;

R = 8.314/m_mol;    %Specific Gas Constant [SPAD, eq 7 .72]
%etac is the combustion efficiency, usually about 0.95 [SPAD]
c_star = etac*sqrt(gamma*R*T_flame)/(gamma*(2/(gamma+1))^((gamma+1)/(2*gamma-2))); %characteristic velocity [SPAD, eq 7.71]

end
