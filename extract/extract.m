% Created by Shreeyam Kacker 2018-02-16
% ICSEDS-EDP
% Extracts data for interpolation

% Gets P vs ox/fuel ratio, vs Cp/Cv, T, molecular weight

%% Housekeeping

clc;
clear;

%% Constants

n = 20;

%% Initialization

% Find pressures

files = ls('*.txt');
m = size(files, 1);

P_cc_vals = zeros(1, m);

for i = 1:length(P_cc_vals)
    P_cc_vals(i) = str2double(files(i, 5:6));
end

% Create matrices (subscript i for interp)

m_mol_data = zeros(n, m);
gamma_data = zeros(n, m);
T_flame_data = zeros(n, m);

%% Main

for i = 1:length(P_cc_vals)
    text = fileread(files(i, :));

    % Molecular weight
    m_mol = findpropepval(text, 'THE MOLECULAR WEIGHT OF THE MIXTURE IS', 5, 4, 2);
    % CP/CV
    cpcv = findpropepval(text, 'CP/CV', 6, 66, 2);
    % T_flame
    t_f = findpropepval(text, 'T(K)', 4, 67, 2);
    
    % Save to matrices
    m_mol_data(:, i) = m_mol';
    gamma_data(:, i) = cpcv';
    T_flame_data(:, i) = t_f';
    
    % Horrible but quick
    if(i == 1)
        % O_F ratio

        o = findpropepval(text, 'NITROUS OXIDE (LIQUID)', 6, 12, 1);
        f = findpropepval(text, 'NYLON 6   POLYAMIDE', 6, 15, 1);

        OF_vals = o./f;
        OF_vals = OF_vals';
    end
end

function values = findpropepval(text, searchstr, length, offset, skip)
    indices = strfind(text, searchstr);
    
    % Discard exhaust pressures
    indices = indices(1:skip:end);
    
    n = size(indices, 2);
    
    values = zeros(1, n);
    
    for i = 1:n
       % Magic number to account for PROPEP's output format
       value_index = indices(i) + strlength(searchstr) + offset;
       values(i) = str2double(text(value_index:value_index + length));
    end
end

