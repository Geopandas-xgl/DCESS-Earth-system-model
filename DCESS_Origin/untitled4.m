 % function [mpr,mdr] = AtmMet_M(t,MH,AT)

% Input : t - time
%         At - atmospheric tracers
% Output: mpr(1) - MH release to atmosphere, 12C. mol/sec 
%         mpr(2) - MH release to atmosphere, 13C. mol/sec 
%         mpr(3) - MH release to one ocean layer, 12C. mol/sec 
%         mpr(4) - MH release to one ocean layer, 13C. mol/sec
%         mdr(1) - methane oxydation in the atmosphere , 12C, mol/sec
%         mdr(2) - methane oxydation in the atmosphere , 13C, mol/sec


global sy rVa mgt R13pdb mdts2 methc13 MRR

t = 0:100:10000;
t = t*sy;
% 
% pCH4o =  0.72e-6;                           %Pre-industrial methane concentration
% pCH4  =  AT(2,1);
% fatm=1.0039;                                %fractionation in atmospheric oxidation
% M     = (pCH4-pCH4o)/pCH4o;
% 
% 
% RCH4o   =  1/(rVa*mdts2*sy);                %PA decay rate for methane, residence time 8.4 yrs
% RCH4    =  RCH4o*(1-0.78*M/(M+11));   
% 
% nlay  = 30;                                 % number of ocean layers receiving methane 
% 
 if (t/sy>=500)
    MRR  =  (mgt.*MHM(3).*(t./sy-MHM(2)).^4.*exp(-MHM(4).*(t./sy-MHM(2))))./sy; 
   
else
     MRR =0;
 end         

 
 
% 
% % mpr(1) =  MH(1)*MRR;
% mpr(2) =  (MH(1)*MRR*(methc13*1e-3+1))*R13pdb;   % methc13 is del 13c of input methane
% mpr(3) =  (1-MH(1))/nlay*MRR;      
% mpr(4) =  (1-MH(1))/nlay*MRR*R13pdb*(methc13*1e-3+1);
% 
% mdr(1) =  pCH4*RCH4;                        % atmospheric methane oxidation to CO2
% mdr(2) =  (pCH4*RCH4*AT(7,1)/pCH4)*fatm;

