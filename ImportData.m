%% Import data from text file.
%
% Script for importing data from the text file
% THERM.txt present in the same folder

%% Initialize variables.
clc;
clear all;
filename = mfilename('fullpath');
filename = strrep(filename, 'ImportData', 'THERM.txt');
startRow = 15;

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%15s%15s%15s%15s%15s%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
textscan(fileID, '%[^\n\r]', startRow-1, 'WhiteSpace', '', 'ReturnOnError', false, 'EndOfLine', '\r\n');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'TextType', 'string', 'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[1,2,3,4,5,6]
    % Converts text in the input cell array to numbers.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        if (mod(row,4)~=1)
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^[-/+]*\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
        end
    end
end

% Generate an intermediate one-dimensional array
[rownum, colnum] = size(raw);
singlearray = cell(1, rownum*colnum);
k = 1;
for i = 1:rownum
    for j = 1:colnum
        singlearray(1,k) = raw(i,j);
          k = k+1;  
    end
end

% Generate the data table consisting of each row for each separate species
rownum = rownum/4;
colnum = colnum*4;
datatable = cell(rownum, colnum);
k = 1;
for i = 1:rownum
    for j = 1:colnum
        datatable(i,j) = singlearray(k);
          k = k+1;  
    end
end

%% Generate the structure containing the name and coefficients for each species
global species;
for i = 1:rownum
    val = datatable(i,1);     species(i).name = val{1};
    val = datatable(i,7);     species(i).a0_high = val{1};
    val = datatable(i,8);     species(i).a1_high = val{1};
    val = datatable(i,9);     species(i).a2_high = val{1};
    val = datatable(i,10);    species(i).a3_high = val{1};
    val = datatable(i,11);    species(i).a4_high = val{1};
    val = datatable(i,13);    species(i).a5_high = val{1};
    val = datatable(i,14);    species(i).a6_high = val{1};
    val = datatable(i,15);    species(i).a0_low = val{1};
    val = datatable(i,16);    species(i).a1_low = val{1};
    val = datatable(i,17);    species(i).a2_low = val{1};
    val = datatable(i,19);    species(i).a3_low = val{1};
    val = datatable(i,20);    species(i).a4_low = val{1};
    val = datatable(i,21);    species(i).a5_low = val{1};
    val = datatable(i,22);    species(i).a6_low = val{1};
end
disp("The name and 14 coefficients for each species are listed in a structure array.");

%% Display the coefficients of CO, H20, CO2 and H2
print_coefficients("CO");
print_coefficients("H2O");
print_coefficients("CO2");
print_coefficients("H2");

%% Calculate the cp, H, S and G of an example element O2
disp("Example calculations for cp, H, S and G for O2 at 300 K:");
O2_heat_capacity_300_K = find_heat_capacity("O2", 300) + " J/K.mol"
O2_enthalpy_300_K = find_enthalpy("O2", 300) + " J/mol"
O2_entropy_300_K = find_entropy("O2", 300) + " J/K.mol"
O2_gibbs_free_energy_300_K = find_gibbs_free_energy("O2", 300) + " J/mol"

%% Display the temperature dependence of Gibbs Free Energy of reaction
print_G_table();

%% Clear temporary variables
clearvars filename startRow formatSpec fileID dataArray ans raw col numericData rawData row regexstr result singlearray numbers invalidThousandsSeparator thousandsRegExp R;


%% ------function for calculating index in the datatable------ %%
%
% ------input : name of species------ %
function [index] = find_index(name)
global species;
rownum = size(species);     rownum = rownum(2);     name = name + " ";
for i = 1:rownum
    if startsWith(species(i).name, name)
        index = i;
        break;
    end
end
end

%% ------function for displaying the coefficients------ %%
%
% ------input : name of species------ %
function [] = print_coefficients(name)
global species;
index = find_index(name);
disp('The coefficients of ' + name + ' are:');
disp(species(index));
end

%% ------function for calculating the Heat Capacity cp------ %%
%
% ------inputs : name of species, temperature------ %
function [cp] = find_heat_capacity(name, T)
global species;
i = find_index(name);
if (T > 1000)
    cp = (species(i).a0_high + (species(i).a1_high*T) + (species(i).a2_high*T*T) + (species(i).a3_high*T*T*T) + (species(i).a4_high*T*T*T*T))*8.31446;
else
    cp = (species(i).a0_low + (species(i).a1_low*T) + (species(i).a2_low*T*T) + (species(i).a3_low*T*T*T) + (species(i).a4_low*T*T*T*T))*8.31446;
end
end

%% ------function for calculating the Enthalpy H------ %%
%
% ------inputs : name of species, temperature------ %
function [H] = find_enthalpy(name, T)
global species;
i = find_index(name);
if (T > 1000)
    H = (species(i).a5_high + (species(i).a0_high*T) + (species(i).a1_high*T*T/2) + (species(i).a2_high*T*T*T/3) + (species(i).a3_high*T*T*T*T/4) + (species(i).a4_high*T*T*T*T*T/5))*8.31446;
else
    H = (species(i).a5_low + (species(i).a0_low*T) + (species(i).a1_low*T*T/2) + (species(i).a2_low*T*T*T/3) + (species(i).a3_low*T*T*T*T/4) + (species(i).a4_low*T*T*T*T*T/5))*8.31446;
end
end

%% ------function for calculating the Entropy S------ %%
%
% ------inputs : name of species, temperature------ %
function [S] = find_entropy(name, T)
global species;
i = find_index(name);
if (T > 1000)
    S = (species(i).a6_high + (species(i).a0_high*log(T)) + (species(i).a1_high*T) + (species(i).a2_high*T*T/2) + (species(i).a3_high*T*T*T/3) + (species(i).a4_high*T*T*T*T/4))*8.31446;
else
    S = (species(i).a6_low + (species(i).a0_low*log(T)) + (species(i).a1_low*T) + (species(i).a2_low*T*T/2) + (species(i).a3_low*T*T*T/3) + (species(i).a4_low*T*T*T*T/4))*8.31446;
end
end

%% ------function for calculating the Gibbs Free Energy G------ %%
%
% ------inputs : name of species, temperature------ %
function [G] = find_gibbs_free_energy(name, T)
H = find_enthalpy(name, T);
S = find_entropy(name, T);
G = H - (T*S);
end

%% ------function for displaying the temperature dependance of Gibbs Free Energy of reaction------ %%
%
% ------input : (no input)------ %
function [] = print_G_table()
global species;
disp("Table showing temperature dependance of the reaction Gibbs Free Energy:");
displaytable = table();
for T = 300:100:1400
    G_reaction = find_gibbs_free_energy("CO2", T) + find_gibbs_free_energy("H2", T) - find_gibbs_free_energy("CO", T) - find_gibbs_free_energy("H2O", T);
    newrow = {T, G_reaction};
    displaytable = [displaytable; newrow];
end
displaytable.Properties.VariableNames = {'Temperature (K)' 'Gibbs Free Energy (J/mol)'};
displaytable
end
