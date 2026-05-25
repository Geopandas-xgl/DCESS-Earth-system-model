function [fanox_global, omz_global]= DCESS_fanox(run)
%% 
% Data ordered accordingly:
% sLL(tracer,level,time) contains low latitude ocean data
% sHL(tracer,level,time) contains high latitude ocean data
% sAT(tracer,zone,time) contains atmospheric data
% st(time) contains model time
% Zones are 
% 1:low latitude (LL), 2:high latitude (HL)
% Ocean tracers are 
% 1:T, 2:S, 3:PO4, 4:DIC, 5:DI13C, 6:DI14C, 7:Alk, 8:O, 9:18O, 
%10:CH4, 11:13CH4, 12: NO3, 13:NH_4, 14:H2S
%Atmosphere tracers are
%1:T, 2:pCH4, 3:pN20, 4:pCO2, 5:pC(13)O2, 6: pC(14)O2, 7:pC(13)H4, 8:pO2
%Land biomasses are
%1:Leaves, 2: Wood, 3:Litter, 4. Soil

global n fdiv aHL aLL 

load ResThilda_M.mat

% ParVal_M
%%
% Establish the atmospheric temperature profile, 2.order legendre pol. in 
% sine of lat, Ta(lat) = PTa(1) + .5*PTa(2) * ( 3*sin(lat)^2 -1 )
% Coefficients are calculated so that the area weighted mean of the profile
% matches the surface mean temperatures in each sector.

t_kyr = run.state.st/1e3;

for ii=1:length(run.state.sAT(1,1,:))
  CTa(1,:) = [ 1 .5*(sin(fdiv)^2-1)                   ];
  CTa(2,:) = [ 1 .5*(sin(fdiv)-sin(fdiv)^3)/(1-sin(fdiv)) ];
  RTa      = [ run.state.sAT(1,1,ii) run.state.sAT(1,2,ii)]';
  PTa(1:2,ii)= CTa\RTa;
end
%% 
GAw = load('GAw_Eo.txt');
GAc = load('GAc_Eo.txt');
%% Calculate time-dependent changes in the anoxic seafloor area based on dissolved oxygen levels
%**********************************************************************************************%
% Setting initial parameters 
OLL = (LL(8,:)*1000)'; % 18O in LL mmol/m3
OHL = (HL(8,:)*1000)'; % 18O in HL mmol/m3
crit = 8.9; % the threshold of O2min for anoxic condition, is equal to that of denitrification (mmol/m3)  
fanox_init = 0.002/(8.9/crit); % The proportion of anoxic seafloor area in the total seafloor area  
fHL = aHL/(aHL+aLL); % the area proportion of high latitude in ocean area
fLL = aLL/(aHL+aLL); % % the area proportion of low latitude in ocean area

% Calculate proportion of each layer in total seafloor area
RS_LL = nan(55,1);
RS_HL =  nan(55,1);
RS_LL(1) = 0;
RS_HL(1) = 0;
for i = 1:54
   RS_LL(i+1) = (GAw(i)-GAw(i+1))/GAw(1); 
   RS_HL(i+1) = (GAc(i)-GAc(i+1))/GAc(1);
end

% Calculate k1 values
syms k1 % variability rate（k < 0）
eqn = fLL*sum(exp(k1.*(OLL/crit-1)).*RS_LL) + fHL*sum(exp(k1.*(OHL/crit-1)).*RS_HL) == fanox_init; % a negative exponential relationship between dissolved oxygen concentration and seafloor anoxia 
k1 = double(vpasolve(eqn, k1));  

% Extract the dissolved oxygen concentration 
ZLL = permute(run.state.sLL(8,:,:),[2,3,1])*1000; % dissolved oxygen concentration in LL mmol/m3
ZHL = permute(run.state.sHL(8,:,:),[2,3,1])*1000; % dissolved oxygen concentration in HL mmol/m3

% Create the matrix of seafloor fraction for each layer
RS_LL_matr = repmat(RS_LL,1,length(t_kyr));
RS_HL_matr = repmat(RS_HL,1,length(t_kyr));

% Transfer dissolved oxygen concentration into anoxic seafloor area 
ZLL_1 = ZLL;
ZLL_1(ZLL_1<=crit) = crit; % The O2min are utilized to repalce the dissolved oxygenation concentration if the dissolved oxygenation concentration are less than O2min     
deg_LL = exp(k1.*(ZLL_1/crit-1)); % Transfer dissolved oxygen concentration into seaflooor anoxic degree for each layer  
ZHL_1 = ZHL;
ZHL_1(ZHL_1<=crit) = crit;
deg_HL = exp(k1.*(ZHL_1/crit-1));
fanox_LL = sum(deg_LL.*RS_LL_matr)*100; % Calculate seafloor aoxic area in LL （unit: %）
fanox_HL = sum(deg_HL.*RS_HL_matr)*100; % Calculate seafloor aoxic area in HL unit: %）
fanox_global = fLL*fanox_LL + fHL*fanox_HL; % % Calculate global seafloor aoxic area in HL unit: %）

%**************************************************************************************************%
%% Calculate time-dependent changes in the OMZ based on dissolved oxygen levels
omz_init = 0.01;
omz_crit = 20;

% Calculate proportion of each layer in total volume
RV_LL = (GAw/GAw(1)*(1/n))';
RV_HL = (GAc/GAc(1)*(1/n))';

% Calculate k1 values
syms k2 % variability rate（k2 < 0）
eqn = fLL*sum(exp(k2.*(OLL/omz_crit-1)).*RV_LL) + fHL*sum(exp(k2.*(OHL/omz_crit-1)).*RV_HL) == omz_init; % a negative exponential relationship between dissolved oxygen concentration and seafloor anoxia 
k2 = double(vpasolve(eqn, k2));

% Create the matrix of area fraction for each layer
RA_LL_matr = repmat(RV_LL,1,length(t_kyr));
RA_HL_matr = repmat(RV_HL,1,length(t_kyr));

% Transfer dissolved oxygen concentration into OMZ
ZLL_2 = ZLL;
ZLL_2(ZLL_2 <= omz_crit) = omz_crit;
deg_LL2 = exp(k2.*(ZLL_2/omz_crit-1)); % Transfer dissolved oxygen concentration into OMZ degree for each layer  
ZHL_2 = ZHL;
ZHL_2(ZHL_2 <= omz_crit) = omz_crit;
deg_HL2 = exp(k2.*(ZHL_2/omz_crit-1));

omz_LL = sum(deg_LL2.*RA_LL_matr)*100; % Calculate seafloor aoxic area in LL （unit: %）
omz_HL = sum(deg_HL2.*RA_HL_matr)*100; % Calculate seafloor aoxic area in HL unit: %）
omz_global = fLL*omz_LL + fHL*omz_HL; % % Calculate global seafloor aoxic area in HL unit: %）
end
