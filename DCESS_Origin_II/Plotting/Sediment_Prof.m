function Sediment_Prof(name)

global sy

sdwpCalA  = ncread(name,'sdwpcala' ); dwpCalA  = squeeze(sdwpCalA(:,:,end));
sdwpCalP  = ncread(name,'sdwpcalp' ); dwpCalP  = squeeze(sdwpCalP(:,:,end));
sdwpCalAr = ncread(name,'sdwpcalar'); dwpCalAr = squeeze(sdwpCalAr(:,:,end));
sdwpCalSO = ncread(name,'sdwpcalso'); dwpCalSO = squeeze(sdwpCalSO(:,:,end));

sdwpOrgA  = ncread(name,'sdwporga' ); dwpOrgA  = squeeze(sdwpOrgA(:,:,end));
sdwpOrgP  = ncread(name,'sdwporgp' ); dwpOrgP  = squeeze(sdwpOrgP(:,:,end));
sdwpOrgAr = ncread(name,'sdwporgar'); dwpOrgAr = squeeze(sdwpOrgAr(:,:,end));
sdwpOrgSO = ncread(name,'sdwporgso'); dwpOrgSO = squeeze(sdwpOrgSO(:,:,end));

swsedA  = ncread(name,'swseda' ); wsedA  = squeeze(swsedA(:,:,end));
swsedP  = ncread(name,'swsedp' ); wsedP  = squeeze(swsedP(:,:,end));
swsedAr = ncread(name,'swsedar'); wsedAr = squeeze(swsedAr(:,:,end));
swsedSO = ncread(name,'swsedso'); wsedSO = squeeze(swsedSO(:,:,end));

ax1    = get_ax;
ax2    = get_ax;

col = '#FF7A28';
Plot_Mod(ax1,dwpCalA,dwpCalP,dwpCalSO,dwpCalAr,1,col,[0 1],'(OM/CaCO_{3})_{dwf}')
col = '#28adff';
Plot_Mod(ax1,dwpOrgA,dwpOrgP,dwpOrgSO,dwpOrgAr,1,col,[0 1],'(OM/CaCO_{3})_{dwf}')

f = sy*1e3;
Plot_Mod(ax2,wsedA,wsedP,wsedSO,wsedAr,f,'k',[0 20],'w_{sed} (cm kyr^{-1})')
end



function Plot_Mod(ax,vA,vP,vSO,vAr,f,col,limits,xlab)
global zm
z   = zm*1e-3;
lw  = 3;
% col = '#386cb0';

msa = vA(5,:);
esa = vA(4,:);
ena = vA(3,:);
mna = vA(2,:);
hna = vA(1,:);

msp = vP(5,:);
esp = vP(4,:);
enp = vP(3,:);
mnp = vP(2,:);
hnp = vP(1,:);

soc = vSO;
arc = vAr;

k=0;

k=k+1;
axes(ax(k))
plot(msa*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(esa*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(ena*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(mna*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(msp*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(esp*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(enp*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(mnp*f,z,'-','Color',col,'LineWidth',lw); hold on


k=k+1;
axes(ax(k))
plot(soc*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(hna*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(hnp*f,z,'-','Color',col,'LineWidth',lw); hold on

k=k+1;
axes(ax(k))
plot(arc*f,z,'-','Color',col,'LineWidth',lw); hold on


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
