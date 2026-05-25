function run = Thilda_M(k_R)

% Integrates the ordinary differential equations in ODE.m using 
% a fourth order Runge Kutta scheme with fixed timestep (h).
tic
clearvars -except  k_R

% Read in initial conditions
%--------------------------------------------
[LL,HL,AT,LB,LLCO3,HLCO3,LLCO3s,HLCO3s,pCalLL,pCalHL,pOrgLL,pOrgHL,dwpCalLL,dwpCalHL,...
 fimLL,fimHL,dwpCalssLL,dwpCalssHL,dwpOrgLL,dwpOrgssLL,dwpOrgHL,dwpOrgssHL,wsedLL,wsedHL ]=InitCnd_M;

global sy n lmdcorg lmdorg dm d rcp lmdcar nto 
global aLL aHL LfLL LfHL ULL UHL GAw GAc NCFLL NCFHL CAFLL CAFHL TSL TSH NPPIL NPPIH workingstate input mixc13 methc13
%% set initial parameter and inject carbon
ParVal_M

input.time_onset = 50*1e3; 

tend = 1.5e5;   % End-time of integration [yr]
t     = 0;

% tend = 2e2;   % End-time of integration [yr]
% t     = (-41*1e3+ input.time_onset)*sy;

input.tend = tend;
dtout  = 100;     % Output interval [yr]
h      = 1/25;    % Timestep [yr], st val 1/25
ko     = 1/h;     % Multiple of timestep at which CO3 and time dependent sediment model are calculated
dt    = h*sy;    % Timestep [sec]
input.k_A = 1;
input.k_R = k_R;

input.meth.total = input.k_R*input.k_A*1800/2;   % Total release of thermogenic methane [GtC], st. val 1000
input.meth.rate_base = 0;
input.meth.frac   = 0.1;   % Fraction of CH4tot that occurs by an amount of time (tfrac) after initial release (at time tCH4)
input.meth.t_start   = -10*1e3 + input.time_onset;    % see above  [yr], st. val. 500 
input.meth.t_end = 70*1e3 + input.time_onset;
input.meth.tfrac  = 5*1e3;    % see above  [yr], st. val. 3500
[input.meth.A,input.meth.B]  = InitRel(input.meth.total,input.meth.frac,input.meth.tfrac);
input.meth.MHM    = [input.meth.t_start input.meth.t_end input.meth.A input.meth.B];

input.mix.total = (1-input.k_R)*(methc13/mixc13)*input.k_A*1800/2;   % Total release of thermogenic methane [GtC], st. val 1000
input.mix.rate_base = 0;
input.mix.frac   = 0.1;   % Fraction of CH4tot that occurs by an amount of time (tfrac) after initial release (at time tCH4)
input.mix.t_start   = -10*1e3 + input.time_onset;    % see above  [yr], st. val. 500 
input.mix.t_end = 70*1e3 + input.time_onset;
input.mix.tfrac  = 5*1e3;    % see above  [yr], st. val. 3500
[input.mix.A,input.mix.B]  = InitRel(input.mix.total,input.mix.frac,input.mix.tfrac);
input.mix.MHM    = [input.mix.t_start input.mix.t_end input.mix.A input.mix.B];

input.mett.total = input.k_A*2700/2;   % Total release of thermogenic methane [GtC], st. val 1000 
input.mett.rate_base = 0;
input.mett.frac   = 0.45;   % Fraction of CH4tot that occurs by an amount of time (tfrac) after initial release (at time tCH4)                         
input.mett.t_start   = -5e3 + input.time_onset;    % see above  [yr], st. val. 500  
input.mett.t_end = 150*1e3 + input.time_onset;
input.mett.tfrac  = 35*1e3;    % see above  [yr], st. val. 3500
[input.mett.A,input.mett.B]  = InitRel(input.mett.total,input.mett.frac,input.mett.tfrac);
input.mett.MHM    = [input.mett.t_start input.mett.t_end input.mett.A input.mett.B];

input.lip.total = input.k_A*1000/2;   % Total release of thermogenic methane [GtC], st. val 1000 
input.lip.rate_base = 0;
input.lip.frac   = 0.1;   % Fraction of CH4tot that occurs by an amount of time (tfrac) after initial release (at time tCH4)                         
input.lip.t_start   = -40*1e3 + input.time_onset;    % see above  [yr], st. val. 500  
input.lip.t_end = 20*1e3 + input.time_onset;
input.lip.tfrac  = 10*1e3;    % see above  [yr], st. val. 3500
[input.lip.A,input.lip.B]  = InitRel(input.lip.total,input.lip.frac,input.lip.tfrac);
input.lip.MHM    = [input.lip.t_start input.lip.t_end input.lip.A input.lip.B];

input.meth.mfrac = 0;
input.mett.mfrac  = 0.9;

filed_names = ["meth","mett","lip","mix"];
t_serial = 0:0.1*1e3:400*1e3; % unit kyr
for i = 1:length(t_serial)
    for j = 1:length(filed_names)
        fname = filed_names(j);
        if ((t_serial(i) >= input.(fname).MHM(1) && t_serial(i) <= input.(fname).MHM(2) ))
             input.(fname).flux(i) = ((input.(fname).MHM (3)*(t_serial(i)-input.(fname).MHM (1))^4*exp(-input.(fname).MHM (4)*(t_serial(i)-input.(fname).MHM (1)))) + input.(fname).rate_base); 
        else
             input.(fname).flux(i) = input.(fname).rate_base;  
        end
    end
end

input.meth.total_serial = cumtrapz(t_serial,input.meth.flux);
input.mett.total_serial = cumtrapz(t_serial,input.mett.flux);
input.lip.total_serial = cumtrapz(t_serial,input.lip.flux);
input.mix.total_serial = cumtrapz(t_serial, input.mix.flux);
%% 
figure("Position",[0,0,500,450])
[ha,~] = tight_subplot(2,1,[0.15,0],[0.15,0.02],[0.1,0.02]);
axes(ha(1))
plot(t_serial,input.meth.flux*2,"LineWidth",1); hold on
plot(t_serial,input.mett.flux*2,"LineWidth",1);
plot(t_serial,input.lip.flux*2,"LineWidth",1);
plot(t_serial,input.mix.flux*2,"LineWidth",1);
ylimits = ylim;
line([input.time_onset,input.time_onset],[ylimits(1),ylimits(2)],"LineStyle","--","Color","k","LineWidth",1);
xlim([0,250e3]);
ticklabels("x",xticks,2,"F");
xtickangle(0);
text_norm(0.7,0.9,"Meth: "+ round(input.meth.total_serial(end)*2,0),10);
text_norm(0.7,0.70,"Mett: "+ round(input.mett.total_serial(end)*2,0),10);
text_norm(0.7,0.50,"Lip: "+ round(input.lip.total_serial(end)*2,0),10);
text_norm(0.7,0.40,"Mix: "+ round(input.mix.total_serial(end)*2,0),10);
xlabel('Time (kyr)');
ylabel("Carbon input (Gtyr^{-1})");
box off
BorderLine(gca,1);
set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.01,0.015],"Layer","top");

axes(ha(2))
data1 = readtable("Data\Xiong-2023-PETM_1.xlsx","Sheet","C and O isotope");
scatter(data1.age_westerhold(1:end-3)*1000+input.time_onset,data1.dC(1:end-3)-mean(data1.dC(1:60)),7,"Marker","o","MarkerEdgeColor","#404040","MarkerFaceColor","#404040","LineWidth",0.01);
xlim([0,250e3]);
xtickangle(0);
xlabel('Time (kyr)');
ylim([-7.5,3]);
yticks(-7:1:3);
ylimits = ylim;
line([input.time_onset,input.time_onset],[ylimits(1),ylimits(2)],"LineStyle","--","Color","k","LineWidth",1);
ticklabels("y",yticks,2,"F");
ylabel('\Delta\delta^{13}C (‰)');
limY = get(gca, 'Ylim');
box off
BorderLine(gca,1);
set(gca,"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.01,0.015],"Layer","top");
print(gcf,"Figure\Input","-dpng","-r600")
%% 
%--------------------------------------------
% Integration
%------------------------------
    for c =1:tend/dtout
  for cc=0:dt:(dtout*sy-dt)
%%---------------
    if mod(t/sy,ko*h)==0
      %---------------
      % Sediment model
      %---------------
    
      % Iterate the carbonate system for the vertical profile of 
      % the CO3 concentration and saturation state 

      [LLCO3,LLCO3s,LLCO2,LLHCO3,LLHp]=CarSysPres_M(LL);
      [HLCO3,HLCO3s,HLCO2,HLHCO3,HLHp]=CarSysPres_M(HL);
      
      [QEBLL,QEBHL,QLL,QHL,aHLNI,Fw,Td] = AtmEnerBal_M(LL(1,1),HL(1,1),AT);
      [RcarLL,RcarHL,RorgLL,RorgHL,Wcarb,Wsil,Worg]=ExtForce_M(AT);
      % Get surface CO2
      [asLL,P14C,CO2LL,CO3LL,HCO3LL,GAMLL] = GasExc_M(LL,AT,ULL,aLL);
      [asHL,P14C,CO2HL,CO3HL,HCO3HL,GAMHL] = GasExc_M(HL,AT,UHL,aHLNI);
      
      % Get New Production
      
      [srcLL,rpl,fnpl,fcpl,fCO2l,fCO3l,fHCO3l] = OrgFlx_M(TSL,NPPIL,LL,AT,aLL  ,LfLL,1,CO2LL,CO3LL,HCO3LL,GAMLL,GAw,pCalLL,pOrgLL,RcarLL,RorgLL,1);
      [srcHL,rph,fnph,fcph,fCO2h,fCO3h,fHCO3h] = OrgFlx_M(TSH,NPPIH,HL,AT,aHLNI,LfHL,0,CO2HL,CO3HL,HCO3HL,GAMHL,GAc,pCalHL,pOrgHL,RcarHL,RorgHL,0);
               
      NPL = (srcLL(3,1)-RorgLL)/aLL;
      NPHM = (srcHL(3,1)-RorgHL)/aHL;
            
      %Time dependent sediment model,
      
      [pCalLL,pOrgLL,dwpCalLL,dwpOrgLL,fimLL,wsedLL]=...
         SMtdCorgNew2_M(TSL,LL,rpl,NCFLL,CAFLL,LLCO3,LLCO3s,NPL,dwpCalLL,dwpCalssLL,dwpOrgLL,dwpOrgssLL,fimLL,wsedLL,ko,h);
      [pCalHL,pOrgHL,dwpCalHL,dwpOrgHL,fimHL,wsedHL]=...
         SMtdCorgNew2_M(TSH,HL,rph,NCFHL,CAFHL,HLCO3,HLCO3s,NPHM,dwpCalHL,dwpCalssHL,dwpOrgHL,dwpOrgssHL,fimHL,wsedHL,ko,h);
      
      RMorg  =  exp(-lmdorg* ( ( dm:d:n ...
          *d )-dm)); 
      RMcorg =  exp(-lmdcorg*( ( dm:d:n*d )-dm));
      RMcar  =  exp(-lmdcar* ( ( dm:d:n*d )-dm));
     
                 
      RMorgLL  = aLL*NPL*       ...
               [-RMorg(1) (RMorg(1:n-2)-RMorg(2:n-1)).*GAw(2:n-1)+...
               RMorg(1:n-2).*(GAw(1:n-2)-GAw(2:n-1)).*pOrgLL(1:n-2) ...
               (RMorg(n-1)-RMorg(n)).*GAw(n)+...
               RMorg(n-1)*(GAw(n-1)-GAw(n))*pOrgLL(n-1)+...
               RMorg(n)*GAw(n)*pOrgLL(n)];
           
      RMorgHL  = aHL*NPHM*       ...
               [-RMorg(1) (RMorg(1:n-2)-RMorg(2:n-1)).*GAc(2:n-1)+...
               RMorg(1:n-2).*(GAc(1:n-2)-GAc(2:n-1)).*pOrgHL(1:n-2) ...
               (RMorg(n-1)-RMorg(n)).*GAc(n)+...
               RMorg(n-1)*(GAc(n-1)-GAc(n))*pOrgHL(n-1)+...
               RMorg(n)*GAc(n)*pOrgHL(n)];
           
      RMcorgLL  = aLL*NPL*rcp*       ...
               [-RMcorg(1) (RMcorg(1:n-2)-RMcorg(2:n-1)).*GAw(2:n-1)+...
               RMcorg(1:n-2).*(GAw(1:n-2)-GAw(2:n-1)).*pOrgLL(1:n-2) ...
               (RMcorg(n-1)-RMcorg(n)).*GAw(n)+...
               RMcorg(n-1)*(GAw(n-1)-GAw(n))*pOrgLL(n-1)+...
               RMcorg(n)*GAw(n)*pOrgLL(n)];
           
      RMcorgHL  = aHL*NPHM*rcp*       ...
               [-RMcorg(1) (RMcorg(1:n-2)-RMcorg(2:n-1)).*GAc(2:n-1)+...
               RMcorg(1:n-2).*(GAc(1:n-2)-GAc(2:n-1)).*pOrgHL(1:n-2) ...
               (RMcorg(n-1)-RMcorg(n)).*GAc(n)+...
               RMcorg(n-1)*(GAc(n-1)-GAc(n))*pOrgHL(n-1)+...
               RMcorg(n)*GAc(n)*pOrgHL(n)];
           
      RMcarLL  = aLL*NPL*rcp*rpl*...
               [-RMcar(1) (RMcar(1:n-2)-RMcar(2:n-1)).*GAw(2:n-1)+...
               RMcar(1:n-2).*(GAw(1:n-2)-GAw(2:n-1)).*pCalLL(1:n-2) ...
               (RMcar(n-1)-RMcar(n)).*GAw(n)+...
               RMcar(n-1)*(GAw(n-1)-GAw(n))*pCalLL(n-1)+...
               RMcar(n)*GAw(n)*pCalLL(n)];
           
      RMcarHL  = aHL*NPHM*rcp*rph*...
               [-RMcar(1) (RMcar(1:n-2)-RMcar(2:n-1)).*GAc(2:n-1)+...
               RMcar(1:n-2).*(GAc(1:n-2)-GAc(2:n-1)).*pCalHL(1:n-2) ...
               (RMcar(n-1)-RMcar(n)).*GAc(n)+...
               RMcar(n-1)*(GAc(n-1)-GAc(n))*pCalHL(n-1)+...
               RMcar(n)*GAc(n)*pCalHL(n)];
                      
      BurorgLL= sum(RMorgLL);
      BurorgHL= sum(RMorgHL);        
      BurcorgLL= sum(RMcorgLL);
      BurcorgHL= sum(RMcorgHL);
      BurcarLL= sum(RMcarLL);
      BurcarHL= sum(RMcarHL); 
           
    end  
   
    t=t+dt;

    [k1LL,k1HL,k1AT,k1LB,srLL1,srHL1,car1] =...
        ODE_M(t,LL,HL,input,AT,LB,pCalLL,pCalHL,pOrgLL,pOrgHL);
      nLL = LL+k1LL*dt/2; 
      nHL = HL+k1HL*dt/2;
      nAT = AT+k1AT*dt/2;
      nLB = LB+k1LB*dt/2;
    [k2LL,k2HL,k2AT,k2LB,srLL2,srHL2,car2] =...
        ODE_M(t+dt/2,nLL,nHL,input,nAT,nLB,pCalLL,pCalHL,pOrgLL,pOrgHL);
      nLL = LL+k2LL*dt/2; 
      nHL = HL+k2HL*dt/2; 
      nAT = AT+k2AT*dt/2;
      nLB = LB+k2LB*dt/2;
    [k3LL,k3HL,k3AT,k3LB,srLL3,srHL3,car3] =...
        ODE_M(t+dt/2,nLL,nHL,input,nAT,nLB,pCalLL,pCalHL,pOrgLL,pOrgHL);
      nLL = LL+k3LL*dt; 
      nHL = HL+k3HL*dt; 
      nAT = AT+k3AT*dt; 
      nLB = LB+k3LB*dt;
    [k4LL,k4HL,k4AT,k4LB,srLL4,srHL4,car4] =...
        ODE_M(t+dt,nLL,nHL,input,nAT,nLB,pCalLL,pCalHL,pOrgLL,pOrgHL);
 

    LL = LL+1/6*dt*( k1LL+2*k2LL+2*k3LL+k4LL );
    HL = HL+1/6*dt*( k1HL+2*k2HL+2*k3HL+k4HL );
    AT = AT+1/6*dt*( k1AT+2*k2AT+2*k3AT+k4AT );
    LB = LB+1/6*dt*( k1LB+2*k2LB+2*k3LB+k4LB );
    
   end % counter cc
  car = 1/6*(car1+2*car2+2*car3+car4);
  meth_vel_12C = car(1);
  mett_vel_12C = car(2);
%% 
  workingstate.scar(:,c) = car; 
  workingstate.st(  c    )   = t/sy;
  workingstate.sAT( :,:,c)   = AT;  
  workingstate.sLL( :,:,c)   = LL;
  workingstate.sHL( :,:,c)   = HL;
  workingstate.sMT( 1,:,c)   = .84*workingstate.sLL(1,:,c)+ .16*workingstate.sHL(1,:,c);

  workingstate.sLB( :,:,c)   = LB; 
  workingstate.srLL(:,:,c)   = 1/6*(srLL1+2*srLL2+2*srLL3+srLL4);
  workingstate.srHL(:,:,c)   = 1/6*(srHL1+2*srHL2+2*srHL3+srHL4);
  
  workingstate.dwcLL(:,c)    = dwpCalLL;
  workingstate.dwcHL(:,c)    = dwpCalHL;
  workingstate.dworgCLL(:,c) = dwpOrgLL;
  workingstate.dworgCHL(:,c) = dwpOrgHL;
  workingstate.fiLL(:,c)     = fimLL;
  workingstate.fiHL(:,c)     = fimHL;
  workingstate.pCalCLL(:,c)  = pCalLL;
  workingstate.pCalCHL(:,c)  = pCalHL;
  workingstate.pOrgCLL(:,c)  = pOrgLL;
  workingstate.pOrgCHL(:,c)  = pOrgHL;
  workingstate.wsLL(:,c)     = wsedLL;
  workingstate.wsHL(:,c)     = wsedHL;
  
  workingstate.QeLL(:,c)     = QEBLL;
  workingstate.QeHL(:,c)     = QEBHL;
  workingstate.ASLL(:,c)     = asLL;
  workingstate.ASHL(:,c)     = asHL;
  
  workingstate.RoLL(1,1,c)   = RorgLL;
  workingstate.RcLL(1,1,c)   = RcarLL;
  workingstate.RoHL(1,1,c)   = RorgHL;
  workingstate.RcHL(1,1,c)   = RcarHL;
  
  workingstate.BurcLL(:,c)   = BurcorgLL;
  workingstate.BurcHL(:,c)   = BurcorgHL;
  workingstate.BurcaLL(:,c)  = BurcarLL;
  workingstate.BurcaHL(:,c)  = BurcarHL;
  workingstate.BurpLL(:,c)   = BurorgLL;
  workingstate.BurpHL(:,c)   = BurorgHL;

  disp(strcat(['Integration at year ' num2str(c*dtout) ...      
		   ' (tend = ',num2str(tend) ')']))
 
end % counter c
%------------------------------
toc
%%
% [fanox_global, omz_global] = DCESS_fanox;
% workingstate.fanox_global = fanox_global;
% workingstate.omz_global = omz_global;
%% 
run.state = workingstate;
run.pars = input;
% save OutThilda_M workingstate 
%% 
% save("dworkspace.mat")
% msgbox('Program has finished running!', 'Notification');
% beep;
end

%------------------------------
%
% Copyright � 2017 Danish Center for Earth System Science
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software 
% and associated documentation files (the "Software"), to deal in the Software without restriction, 
% including without limitation the rights to use, copy, and modify copies of the Software, and 
% to permit persons to whom the Software is furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all copies or 
% substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING 
% BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, 
% DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, 
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

