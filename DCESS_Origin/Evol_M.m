function Evol_M

close all
clear all

% Plot tracer evolutions
%-----
load('OutThilda_M')

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

%-----

figure(1); clf

set(gcf,'paperunits','inches','units','inches')
set(gcf,'paperposition',[1 2 6 7])
set(gcf,'position',     [0 0 6 7])
set(0  ,'defaultlinelinewidth',2)
set(0  ,'defaultaxesfontname','times')
set(0  ,'defaultaxesfontsize',8)
%-----

% Read in parameter values
%-----
ParVal_M
global sy rcp n fdiv R13pdb rpm

for ii=2:length(sAT(4,1,:))
    st(ii)=100*ii;
           
    [pCO2,K0,CO2,CO3,HCO3,GAM]=CarSys_M(sLL(1,1,ii),sLL(2,1,ii),sLL(4,1,ii),sLL(7,1,ii));
    GAMLL(1,1,ii)=GAM;
    [pCO2,K0,CO2,CO3,HCO3,GAM]=CarSys_M(sHL(1,1,ii),sHL(2,1,ii),sHL(4,1,ii),sHL(7,1,ii));
    GAMHL(1,1,ii)=GAM;
     
       for i=2:n
         if dwcLL(i,ii)<0.1
             CCDLL10(ii)=(i-1)*100+(dwcLL(i-1,ii)-0.1)./(dwcLL(i-1,ii)-dwcLL(i,ii))*100; break
         end
        end
       for i=2:n
         if dwcHL(i,ii)<0.1
             CCDHL10(ii)=(i-1)*100+(dwcHL(i-1,ii)-0.1)./(dwcHL(i-1,ii)-dwcHL(i,ii))*100; break
         else
             CCDHL10(ii)=5500;
         end
        end
end     

% Calculate delta13 and Delta 14 values in permil for ocean and atmosphere
%-----

d13a   = squeeze( (sAT(5,1,:)./sAT(4,1,:)/R13pdb-1)*1e3 );
d13ao   = squeeze( (sAT(5,1,1)./sAT(4,1,1)/R13pdb-1)*1e3 );      
d13ad = d13a - d13ao;                                            %atmosphere

d13LL100mo = squeeze( (sLL(5,1,1)./sLL(4,1,1)/R13pdb-1)*1e3 );
d13LL1000mo = squeeze( (sLL(5,10,1)./sLL(4,10,1)/R13pdb-1)*1e3 );
d13LL3000mo = squeeze( (sLL(5,30,1)./sLL(4,30,1)/R13pdb-1)*1e3 );

d13LL100m = squeeze( (sLL(5,1,:)./sLL(4,1,:)/R13pdb-1)*1e3 );
d13LL1000m = squeeze( (sLL(5,10,:)./sLL(4,10,:)/R13pdb-1)*1e3 );
d13LL3000m = squeeze( (sLL(5,30,:)./sLL(4,30,:)/R13pdb-1)*1e3 );

d13LL100md = d13LL100m - d13LL100mo;                       %ocean, 100m
d13LL1000md = d13LL1000m - d13LL1000mo;                    %ocean, 1000m
d13LL3000md = d13LL3000m - d13LL3000mo;                    %ocean, 3000m 

d18cLL100mo = squeeze( sLL(9,1,1).*1000+(16.5-sLL(1,1,1))./4.8 );
d18cLL1000mo = squeeze( sLL(9,10,1).*1000+(16.5-sLL(1,10,1))./4.8 );
d18cLL3000mo = squeeze( sLL(9,30,1).*1000+(16.5-sLL(1,30,1))./4.8 );

d18cLL100m = squeeze( sLL(9,1,:).*1000+(16.5-sLL(1,1,:))./4.8 );
d18cLL1000m = squeeze( sLL(9,10,:).*1000+(16.5-sLL(1,10,:))./4.8 );
d18cLL3000m = squeeze( sLL(9,30,:).*1000+(16.5-sLL(1,30,:))./4.8 );

d18cLL100md = d18cLL100m-d18cLL100mo;                      %ocean, 100m
d18cLL1000md = d18cLL1000m-d18cLL1000mo;                   %ocean, 1000m
d18cLL3000md = d18cLL3000m-d18cLL3000mo;                   %ocean, 3000m

%-----

subplot(6,2,1)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
% Establish the atmospheric temperature profile, 2.order legendre pol. in 
% sine of lat, Ta(lat) = PTa(1) + .5*PTa(2) * ( 3*sin(lat)^2 -1 )
% Coefficients are calculated so that the area weighted mean of the profile
% matches the surface mean temperatures in each sector.
for ii=1:length(sAT(1,1,:))
  CTa(1,:) = [ 1 .5*(sin(fdiv)^2-1)                   ];
  CTa(2,:) = [ 1 .5*(sin(fdiv)-sin(fdiv)^3)/(1-sin(fdiv)) ];
  RTa      = [ sAT(1,1,ii) sAT(1,2,ii)]';
  PTa(1:2,ii)= CTa\RTa;
end
plot(st/1e3,squeeze(sAT(1,1,:)),'r');hold on
plot(st/1e3,squeeze(sAT(1,2,:)),'b');
plot(st/1e3,PTa(1,:),'k--');
axis([[min(st) max(st)]/1e3 -10 45])
ylabel('T^{atm} (^oC)')

%-------------

subplot(6,2,2)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st/1e3,squeeze(sHL(1,1,:)),'b-'); hold on
plot(st/1e3,squeeze(sLL(1,1,:)),'r-')
plot(st/1e3,0.84*squeeze(sLL(1,1,:))+0.16*squeeze(sHL(1,1,:)),'g-')
ylabel('T^{ocean}(^oC)')
%-------------

subplot(6,2,3)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])

plot(st/1e3,squeeze(sAT(4,1,:))*1e6,'b-'); hold on
ylabel('pCO_2 (ppm)')
%-------------

subplot(6,2,4)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st/1e3,squeeze(sAT(2,1,:))*1e6,'b-'); hold on 
plot(st/1e3,squeeze(sAT(3,1,:))*1e7,'r-')    %N20 x 10 (ppm)for visualization
ylabel('pCH4, pN2O (ppm)')
%-------------

subplot(6,2,5)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])

Tr  = 10;                                 % [oC] Reference temperatur, st. val 10
p1  = 1;                                  % [-] O.Marchal 98,Maier-Reimer 93
p2  = 0.1565;                             % [oC^-1] Maier-Reimer 93; revised 22/08/2014
rpLL(1,1,:)  = rpm* p1*exp(p2*(sLL(1,1,:)-Tr)) ...
    ./(1+p1*exp(p2*(sLL(1,1,:)-Tr))).*(GAMLL(1,1,:)-1)./(1+(GAMLL(1,1,:)-1));  
rpHL(1,1,:)  = rpm* p1*exp(p2*(sHL(1,1,:)-Tr)) ...
    ./(1+p1*exp(p2*(sHL(1,1,:)-Tr))).*(GAMHL(1,1,:)-1)./(1+(GAMHL(1,1,:)-1));
NPLL(1,1,:)=-2*(srLL(3,1,:)-RoLL(1,1,:))*rcp*12.011*sy/1e15;
NPHL(1,1,:)=-2*(srHL(3,1,:)-RoHL(1,1,:))*rcp*12.011*sy/1e15;
plot(st/1e3,squeeze(NPLL(1,1,:).*rpLL(1,1,:)),'r:'); hold on     %LL CarbC
plot(st/1e3,squeeze(NPHL(1,1,:).*rpHL(1,1,:)),'b:')              %HL CarbC 
plot(st/1e3,squeeze(NPLL(1,1,:)),'r')                            %LL OrgC 
plot(st/1e3,squeeze(NPHL(1,1,:)),'b')                            %HL OrgC 
axis([[min(st) max(st)]/1e3 0 6])
ylabel('NP (Gton C/yr)')
%-------------

subplot(6,2,6)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])

plot(st/1e3,squeeze(sLL(8,30,:)),'b-'); hold on
plot(st/1e3,squeeze(sLL(8,10,:)),'g-');
plot(st/1e3,squeeze(sLL(8,1,:)),'r-');
axis([[min(st) max(st)]/1e3 0 0.3])
ylabel('O_2 (mol m^{-3})')
%-------------

subplot(6,2,7)
%-------------
%Excursions are plotted
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st/1e3,d18cLL100md,'r-'); hold on
plot(st/1e3,d18cLL1000md,'g-')
plot(st/1e3,d18cLL3000md,'b-')
ylabel('\delta^{18}O_{carb}(o/oo)')

%-------------

subplot(6,2,8)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st/1e3,2*squeeze(sLB(1,1,:) + sLB(2,1,:)),'b-'); hold on     %above ground
plot(st/1e3,2*squeeze(sLB(3,1,:) + sLB(4,1,:)),'r-');             %below ground
plot(st/1e3,2*squeeze(sLB(1,1,:) + sLB(2,1,:) + sLB(3,1,:) + sLB(4,1,:)),'k-');
ylabel('LB (GtC)')
%-------------

subplot(6,2,9)
%-------------
% Excursions are plotted
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st/1e3,d13ad,'m-');hold on
plot(st/1e3,d13LL100md,'r-'); hold on
plot(st/1e3,d13LL1000md,'g-')
plot(st/1e3,d13LL3000md,'b-')
ylabel('\delta^{13}C_{atm,oc}(o/oo)')
xlabel('t (kyr)')

%-------------
subplot(6,2,10)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st(2:end)/1e3,-CCDLL10(2:end),'r'); hold on
plot(st(2:end)/1e3,-CCDHL10(2:end),'b')
axis([[min(st) max(st)]/1e3 -5500 -1000])
xlabel('t (kyr)')
ylabel('CCD (m)')
%-------------

subplot(6,2,11)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st/1e3,2*(BurcLL+BurcHL),'r'); hold on
plot(st/1e3,2*(BurcaLL+BurcaHL),'b');
ylabel('Org/CalC burial (molC/s)')
xlabel('t (kyr)')

subplot(6,2,12)
%-------------
P=get(gca,'position');set(gca,'position',[P(1) P(2) P(3) P(4)*1.2])
plot(st/1e3,wsLL(10,:)*1e3*sy,'r-'); hold on
plot(st/1e3,wsLL(30,:)*1e3*sy,'r:')
plot(st/1e3,wsHL(10,:)*1e3*sy,'b')
plot(st/1e3,wsHL(30,:)*1e3*sy,'b:')
ylabel('Sed. vel.(cm/kyr)')
xlabel('t (kyr)')

return
