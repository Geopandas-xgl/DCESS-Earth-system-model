function Vertical_Prof(name)


sHLNA = ncread(name,'shlna'); HLNA = squeeze(sHLNA(:,:,end));
sMLNA = ncread(name,'smlna'); MLNA = squeeze(sMLNA(:,:,end));
sLLNA = ncread(name,'sllna'); LLNA = squeeze(sLLNA(:,:,end));
sLLSA = ncread(name,'sllsa'); LLSA = squeeze(sLLSA(:,:,end));
sMLSA = ncread(name,'smlsa'); MLSA = squeeze(sMLSA(:,:,end));

sHLNP = ncread(name,'shlnp'); HLNP = squeeze(sHLNP(:,:,end));
sMLNP = ncread(name,'smlnp'); MLNP = squeeze(sMLNP(:,:,end));
sLLNP = ncread(name,'sllnp'); LLNP = squeeze(sLLNP(:,:,end));
sLLSP = ncread(name,'sllsp'); LLSP = squeeze(sLLSP(:,:,end));
sMLSP = ncread(name,'smlsp'); MLSP = squeeze(sMLSP(:,:,end));

sArc  = ncread(name,'sarc'); Arc = squeeze(sArc(:,:,end));
sSOc  = ncread(name,'ssoc'); SOc = squeeze(sSOc(:,:,end));



% name_file = 'E:\Trabajo\DCESS\DCESS_II\Datos_Obs\GLODAPv2_2022\Analysis_Dec2023\GLODAP_v2022_Dec2023_ModelUnits.mat';
M   = load('GLODAP_v2022_Dec2023_ModelUnits.mat');
O18 = load('O18_Schmidt_v1_22.mat');

% T --> °C
% S --> g kg-1
% PO4, DIC, ALK and O2 --> g kg-1
% δ13C and ∆14C --> permil


% Plot observed data

axT    = get_ax;
axS    = get_ax;
axO18  = get_ax;
axPO4  = get_ax;
axO2   = get_ax;
axDIC  = get_ax;
axALK  = get_ax;
axd13C = get_ax;
axD14C = get_ax;

fvar = 1028*1e-6*1e3; % m mol m-3
Plot_Obs(axT  ,M.MAt.CT ,M.MPa.CT ,M.MSO.CT ,M.MAr.CT ,1   ,[-5   25],'T (°C)');
Plot_Obs(axS  ,M.MAt.SA ,M.MPa.SA ,M.MSO.SA ,M.MAr.SA ,1   ,[33   37],'S (g kg^{-1})');
Plot_Obs_O18(axO18,O18.At,O18.Pa,O18.SO,O18.Ar,[-2 2],'\delta^{18}O_{w} (permil)')
Plot_Obs(axPO4,M.MAt.PO4,M.MPa.PO4,M.MSO.PO4,M.MAr.PO4,fvar,[0    4 ],'PO_{4} (mmol m^{-3})');
Plot_Obs(axO2 ,M.MAt.O2 ,M.MPa.O2 ,M.MSO.O2 ,M.MAr.O2 ,fvar,[0   500],'O_{2} (mmol m^{-3})');
fvar = 1028*1e-6; % m mol m-3
Plot_Obs(axDIC ,M.MAt.DIC ,M.MPa.DIC ,M.MSO.DIC ,M.MAr.DIC ,fvar,[1.9  2.5],'DIC (mol m^{-3})');
Plot_Obs(axALK ,M.MAt.ALK ,M.MPa.ALK ,M.MSO.ALK ,M.MAr.ALK ,fvar,[2.1  2.6],'ALK (mol m^{-3})');
Plot_Obs(axd13C,M.MAt.d13C,M.MPa.d13C,M.MSO.d13C,M.MAr.d13C,1   ,[-1.5 2.5],'\delta^{13}C (permil)');
Plot_Obs(axD14C,M.MAt.D14C,M.MPa.D14C,M.MSO.D14C,M.MAr.D14C,1   ,[-275 0  ],'\Delta^{14}C (permil)');

iT=1;
Plot_Mod(axT,MLSA(iT,:),LLSA(iT,:),LLNA(iT,:),MLNA(iT,:),HLNA(iT,:),...
             MLSP(iT,:),LLSP(iT,:),LLNP(iT,:),MLNP(iT,:),HLNP(iT,:),SOc(iT,:),Arc(iT,:),1);
         
iT=2;
Plot_Mod(axS,MLSA(iT,:),LLSA(iT,:),LLNA(iT,:),MLNA(iT,:),HLNA(iT,:),...
             MLSP(iT,:),LLSP(iT,:),LLNP(iT,:),MLNP(iT,:),HLNP(iT,:),SOc(iT,:),Arc(iT,:),1);

iT=10;
Plot_Mod(axO18,MLSA(iT,:),LLSA(iT,:),LLNA(iT,:),MLNA(iT,:),HLNA(iT,:),...
               MLSP(iT,:),LLSP(iT,:),LLNP(iT,:),MLNP(iT,:),HLNP(iT,:),SOc(iT,:),Arc(iT,:),1e3);

iT=3;
Plot_Mod(axPO4,MLSA(iT,:),LLSA(iT,:),LLNA(iT,:),MLNA(iT,:),HLNA(iT,:),...
               MLSP(iT,:),LLSP(iT,:),LLNP(iT,:),MLNP(iT,:),HLNP(iT,:),SOc(iT,:),Arc(iT,:),1e3);

iT=9;
Plot_Mod(axO2,MLSA(iT,:),LLSA(iT,:),LLNA(iT,:),MLNA(iT,:),HLNA(iT,:),...
              MLSP(iT,:),LLSP(iT,:),LLNP(iT,:),MLNP(iT,:),HLNP(iT,:),SOc(iT,:),Arc(iT,:),1e3);

iT=5;
Plot_Mod(axDIC,MLSA(iT,:),LLSA(iT,:),LLNA(iT,:),MLNA(iT,:),HLNA(iT,:),...
               MLSP(iT,:),LLSP(iT,:),LLNP(iT,:),MLNP(iT,:),HLNP(iT,:),SOc(iT,:),Arc(iT,:),1);
           
iT=8;
Plot_Mod(axALK,MLSA(iT,:),LLSA(iT,:),LLNA(iT,:),MLNA(iT,:),HLNA(iT,:),...
               MLSP(iT,:),LLSP(iT,:),LLNP(iT,:),MLNP(iT,:),HLNP(iT,:),SOc(iT,:),Arc(iT,:),1);

% d13C, D14C
[d13CmsA,D14CmsA]=get_dC13_D14C(MLSA);
[d13CesA,D14CesA]=get_dC13_D14C(LLSA);
[d13CenA,D14CenA]=get_dC13_D14C(LLNA);
[d13CmnA,D14CmnA]=get_dC13_D14C(MLNA);
[d13ChnA,D14ChnA]=get_dC13_D14C(HLNA);

[d13CmsP,D14CmsP]=get_dC13_D14C(MLSP);
[d13CesP,D14CesP]=get_dC13_D14C(LLSP);
[d13CenP,D14CenP]=get_dC13_D14C(LLNP);
[d13CmnP,D14CmnP]=get_dC13_D14C(MLNP);
[d13ChnP,D14ChnP]=get_dC13_D14C(HLNP);

[d13CSOc,D14CSOc]=get_dC13_D14C(SOc);
[d13CArc,D14CArc]=get_dC13_D14C(Arc);

Plot_Mod(axd13C,d13CmsA,d13CesA,d13CenA,d13CmnA,d13ChnA,d13CmsP,d13CesP,d13CenP,d13CmnP,d13ChnP,d13CSOc,d13CArc,1);
Plot_Mod(axD14C,D14CmsA,D14CesA,D14CenA,D14CmnA,D14ChnA,D14CmsP,D14CesP,D14CenP,D14CmnP,D14ChnP,D14CSOc,D14CArc,1);

end


% External
function Plot_Obs(ax,vA,vP,vSO,vAr,fvar,limits,xlab)
mkt = 'o';
mks = 5;
mkc = '#5ac0ff';
fz = 1e-3;

k=0;
if strcmp(xlab,'\Delta^{14}C (permil)')
    k=k+1;
    axes(ax(k))
    C14x = vA.ms.var(vA.ms.var<=0)*fvar; zx = vA.ms.z(vA.ms.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vA.es.var(vA.es.var<=0)*fvar; zx = vA.es.z(vA.es.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vA.en.var(vA.en.var<=0)*fvar; zx = vA.en.z(vA.en.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vA.mn.var(vA.mn.var<=0)*fvar; zx = vA.mn.z(vA.mn.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on


    k=k+1;
    axes(ax(k))
    C14x = vP.ms.var(vP.ms.var<=0)*fvar; zx = vP.ms.z(vP.ms.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vP.es.var(vP.es.var<=0)*fvar; zx = vP.es.z(vP.es.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vP.en.var(vP.en.var<=0)*fvar; zx = vP.en.z(vP.en.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vP.mn.var(vP.mn.var<=0)*fvar; zx = vP.mn.z(vP.mn.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on
    
    
    k=k+1;
    axes(ax(k))
    C14x = vSO.var(vSO.var<=0)*fvar; zx = vSO.z(vSO.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vA.hn.var(vA.hn.var<=0)*fvar; zx = vA.hn.z(vA.hn.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vP.hn.var(vP.hn.var<=0)*fvar; zx = vP.hn.z(vP.hn.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    C14x = vAr.var(vAr.var<=0)*fvar; zx = vAr.z(vAr.var<=0)*fz;
    C14 = C14x(zx>=1); z = zx(zx>=1);
    plot(C14,z,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on
else 
    k=k+1;
    axes(ax(k))
    plot(vA.ms.var*fvar,vA.ms.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vA.es.var*fvar,vA.es.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vA.en.var*fvar,vA.en.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vA.mn.var*fvar,vA.mn.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on


    k=k+1;
    axes(ax(k))
    plot(vP.ms.var*fvar,vP.ms.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vP.es.var*fvar,vP.es.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vP.en.var*fvar,vP.en.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vP.mn.var*fvar,vP.mn.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vSO.var*fvar,vSO.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vA.hn.var*fvar,vA.hn.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vP.hn.var*fvar,vP.hn.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

    k=k+1;
    axes(ax(k))
    plot(vAr.var*fvar,vAr.z*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on
end


set(ax(1:12),'Xlim',limits,'Ylim',[0 5.5],'YDir','Reverse')
set(ax(1:4) ,'XAxisLocation','Top');
set(ax(9:12),'XAxisLocation','Bottom')

set(ax([2:4 6:8 10:12]),'YTicklabel',{''})
set(ax(5:8),'XTicklabel',{''})

ylabel(ax([1 5 9]),'Depth (km)')
xlabel(ax(1:4),xlab)
xlabel(ax(9:12),xlab)


txt = {'msA','esA','enA','mnA',...
       'msP','esP','enP','mnP',...
       'SOc','hnA','hnP','Arc'};
for i=1:12
    text(ax(i),.01,.075,txt{i},'Units','normalized','FontWeight','bold','FontSize',12);
end

end

function Plot_Obs_O18(ax,O18A,O18P,O18SO,O18Ar,limits,xlab)

mkt = 'o';
mks = 5;
mkc = '#5ac0ff';
fz = 1e-3;

k=0;
k=k+1;
axes(ax(k))
plot(O18A.ms(:,4),O18A.ms(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18A.es(:,4),O18A.es(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18A.en(:,4),O18A.en(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18A.mn(:,4),O18A.mn(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on


k=k+1;
axes(ax(k))
plot(O18P.ms(:,4),O18P.ms(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18P.es(:,4),O18P.es(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18P.en(:,4),O18P.en(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18P.mn(:,4),O18P.mn(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on


k=k+1;
axes(ax(k))
plot(O18SO(:,4),O18SO(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18A.hn(:,4),O18A.hn(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18P.hn(:,4),O18P.hn(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

k=k+1;
axes(ax(k))
plot(O18Ar(:,4),O18Ar(:,3)*fz,mkt,'MarkerFaceColor',mkc,'MarkerEdgeColor',mkc,'MarkerSize',mks); hold on

set(ax(1:12),'Xlim',limits,'Ylim',[0 5.5],'YDir','Reverse')
set(ax(1:4) ,'XAxisLocation','Top');
set(ax(9:12),'XAxisLocation','Bottom')

set(ax([2:4 6:8 10:12]),'YTicklabel',{''})
set(ax(5:8),'XTicklabel',{''})

ylabel(ax([1 5 9]),'Depth (km)')
xlabel(ax(1:4),xlab)
xlabel(ax(9:12),xlab)


txt = {'msA','esA','enA','mnA',...
       'msP','esP','enP','mnP',...
       'SOc','hnA','hnP','Arc'};
for i=1:12
    text(ax(i),.01,.075,txt{i},'Units','normalized','FontWeight','bold','FontSize',12);
end
end

function Plot_Mod(ax,msa,esa,ena,mna,hna,msp,esp,enp,mnp,hnp,soc,arc,f)
global zm
z   = zm*1e-3;
lw  = 3;
col = '#386cb0';

hna(37:end) = NaN;

k=0;

k=k+1;
axes(ax(k))
plot(msa*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(esa*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(ena*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(mna*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(msp*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(esp*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(enp*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(mnp*f,z,'-','Color',col,'LineWidth',lw)


k=k+1;
axes(ax(k))
plot(soc*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(hna*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(hnp*f,z,'-','Color',col,'LineWidth',lw)

k=k+1;
axes(ax(k))
plot(arc*f,z,'-','Color',col,'LineWidth',lw)
end


function [d13C,D14C]=get_dC13_D14C(Psi)
global R13pdb R14oas

C12 = Psi(5,:);
C13 = Psi(6,:);
C14 = Psi(7,:);

d13C = ( C13./C12/R13pdb - 1 )*1e3;
D14C = ( C14./C12/R14oas.*( R13pdb*0.975./(C13./C12) ).^2 - 1 )*1e3 ;

end























