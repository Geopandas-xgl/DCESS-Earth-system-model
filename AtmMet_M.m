function [mpr,mdr,MRR] = AtmMet_M(t,input ...
    ,AT)

% Input : t - time
%         At - atmospheric tracers
% Output: mpr(1) - MH release to atmosphere, 12C. mol/sec 
%         mpr(2) - MH release to atmosphere, 13C. mol/sec 
%         mpr(3) - MH release to one ocean layer, 12C. mol/sec 
%         mpr(4) - MH release to one ocean layer, 13C. mol/sec
%         mdr(1) - methane oxydation in the atmosphere , 12C, mol/sec
%         mdr(2) - methane oxydation in the atmosphere , 13C, mol/sec


global sy rVa mgt R13pdb mdts2 methc13 mettc13

pCH4o =  0.72e-6;                           %Pre-industrial methane concentration
pCH4  =  AT(2,1);
fatm=1.0039;                                % C fractionation in atmospheric oxidation
M     = (pCH4-pCH4o)/pCH4o;

                                             
RCH4o   =  1/(rVa*mdts2*sy);                %PA decay rate for methane, residence time 8.4 yrs
RCH4    =  RCH4o*(1-0.78*M/(M+11));         % Decay rate is the inverse of residence time 

%************************************
nlay  = 20-9+1;                                 % numbers of ocean layers receiving methane (from 900 to 2000); nlayer = 16; See ODE_M.m 
%*************************************

%------------------------------------------------

%%%%%%% Hydrate C release rate derive from oxidation mol/sy
if ((t/sy >= input.meth.MHM(1) && t/sy <= input.meth.MHM(2) ))
     MRRh = ((input.meth.MHM (3)*(t/sy-input.meth.MHM (1))^4*exp(-input.meth.MHM (4)*(t/sy-input.meth.MHM (1)))) + input.meth.rate_base)*mgt/sy; 
else
     MRRh = input.meth.rate_base*mgt/sy;  
end

if ((t/sy >= input.mett.MHM(1) && t/sy <= input.mett.MHM(2) ))
     MRRt = ((input.mett.MHM (3)*(t/sy-input.mett.MHM (1))^4*exp(-input.mett.MHM (4)*(t/sy-input.mett.MHM (1)))) + input.mett.rate_base)*mgt/sy; 
else
     MRRt = input.mett.rate_base*mgt/sy;  
end


MRR = [MRRh, MRRt];

mpr(1) =  input.meth.mfrac*MRRh + input.mett.mfrac*MRRt;
mpr(2) =  (input.meth.mfrac*MRRh*(methc13*1e-3+1))*R13pdb + (input.mett.mfrac*MRRt*(mettc13*1e-3+1))*R13pdb;   % methc13 is del 13c of input methane
mpr(3) =  (1-input.meth.mfrac)/nlay*MRRh + (1-input.mett.mfrac)/nlay*MRRt;      
mpr(4) =  (1-input.meth.mfrac)/nlay*MRRh*R13pdb*(methc13*1e-3+1)+ (1-input.mett.mfrac)/nlay*MRRt*R13pdb*(mettc13*1e-3+1);

mdr(1) =  pCH4*RCH4;  % atmospheric methane oxidation to CO2; pCH4 is the concentration of CH4 in atmosphere
mdr(2) =  (pCH4*RCH4*AT(7,1)/pCH4)*fatm;

return


