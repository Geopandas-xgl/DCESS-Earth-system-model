function Atm_results(name)

global dr

sAT = ncread(name,'sat');
AT  = squeeze(sAT(:,:,end));


[TamN,TamS,TamG]    = Atm_mean(AT(1,:));
[~,~,pCH4G]   = Atm_mean(AT(3,:));
[~,~,pN2OG]   = Atm_mean(AT(4,:));
[~,~,pCO2G]   = Atm_mean(AT(5,:));
[~,~,p13CO2G] = Atm_mean(AT(6,:));

global R13pdb

d13Ca =  ( p13CO2G./pCO2G/R13pdb - 1 )*1e3;

disp('Global mean Ta [SH NH Global] (°C)')
disp(num2str([TamS TamN TamG]))

disp('Global mean pCH4 (ppm)')
disp(num2str(pCH4G*1e6))

disp('Global mean pN2O (ppm)')
disp(num2str(pN2OG*1e6))

disp('Global mean pCO2 (ppm)')
disp(num2str(pCO2G*1e6))

disp('Global mean d13C_atm (permil)')
disp(num2str(d13Ca))

disp('Sea-ice line position [North Atl North Pac SO] (degrees)')
disp(num2str(AT(10,[1 2 6])/dr))



end




function [vmN,vmS,vmG] = Atm_mean(var)


global fdh fdl fde
vmN = var(1)*(1-sin(fdh)) + var(2)*(sin(fdh)-sin(fdl)) + var(3)*(sin(fdl)-sin(fde)) ;
vmS = var(6)*(1-sin(fdh)) + var(5)*(sin(fdh)-sin(fdl)) + var(4)*(sin(fdl)-sin(fde)) ;
vmG  = .5*(vmN+vmS);

end