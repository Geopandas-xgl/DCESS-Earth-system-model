% Simple file to reproduce manuscript's figures
clear
clc
close all

Globals;

name='../Out.nc';

Atm_results(name);
Oce_results(name);
Lithos_results;

% Ocean circulation
Ocean_Circulation(name)

% Vertical profiles
Vertical_Prof(name);

% Sediment Profiles
Sediment_Prof(name);

