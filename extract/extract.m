% Created by Shreeyam Kacker 2018-02-16 1759
% ICSEDS-EDP
% Extracts data for interpolation

% Gets P vs ox/fuel ratio, vs Cp/Cv, T, molecular weight
% DOES NOT WORK ON MAC OS DUE TO FILE SYSTEM DIFFERENCES

%% Housekeeping

clc;
clear;

%% Constants

n = 20;

%% Initialization

% Get files list

files = ls('*.txt');
m = size(files, 1);

P_cc_vals = zeros(1, m);    % [Pa]

% Find pressures from file inputs
for i = 1:length(P_cc_vals)
    P_cc_vals(i) = str2double(files(i, 5:6)) * 1e5;
end

% Create matrices

m_mol_data = zeros(n, m);   % [kg/mol]
gamma_data = zeros(n, m);   % [-]
T_flame_data = zeros(n, m); % [K]

%% Main

for i = 1:length(P_cc_vals)
    % Read text from file
    text = fileread(files(i, :));

    % Molecular weight
    m_mol = findpropepval(text, 'THE MOLECULAR WEIGHT OF THE MIXTURE IS', 5, 4, 2);
    % CP/CV
    cpcv = findpropepval(text, 'CP/CV', 6, 66, 2);
    % T_flame
    t_f = findpropepval(text, 'T(K)', 4, 67, 2);
    
    % Save to matrices
    m_mol_data(:, i) = m_mol' * 1e-3;  
    gamma_data(:, i) = cpcv';
    T_flame_data(:, i) = t_f';
    
    % Horrible but quick
    if(i == 1)
        % O_F ratio

        o = findpropepval(text, 'NITROUS OXIDE (LIQUID)', 6, 12, 1);
        f = findpropepval(text, 'NYLON 6   POLYAMIDE', 6, 15, 1);

        OF_vals = o./f;
        OF_vals = OF_vals'; % [-]
    end
end

% Finds a numeric values that is offset some distance from the search
% string

% ARGUMENTS:

% Text          Text to search
% Searchstr     Searchstr to look for
% Length        Characterwise width of numeric value
% Offset        Distance from search string to value
% Skip          Discard every (skip)th value

% OUTPUTS:

% Numerical values (double)

function values = findpropepval(text, searchstr, length, offset, skip)
    % Find indices of search strings
    indices = strfind(text, searchstr);
    
    % Discard exhaust pressures, etc
    indices = indices(1:skip:end);
    
    n = size(indices, 2);
    
    values = zeros(1, n);
    
    for i = 1:n
       % Accounting for PROPEP's file format
       value_index = indices(i) + strlength(searchstr) + offset;
       % Save to output values
       values(i) = str2double(text(value_index:value_index + length));
    end
end

