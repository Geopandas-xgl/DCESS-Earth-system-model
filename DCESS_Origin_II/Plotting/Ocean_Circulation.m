function Ocean_Circulation(name)

global n dr zm Bsh Bs fdh fdl fde fBr Bn

qA    = ncread(name,'q_atl');
qP    = ncread(name,'q_pac');
Ekman = ncread(name,'ekman');

EVFsA = Ekman(:,1);
EVFsP = Ekman(:,2);

qA    = [qA(:,1) qA(:,1)+EVFsA qA(:,2) qA(:,3) qA(:,4) qA(:,5)  qA(:,6)]*1e-6;
qP    = [qP(:,1) qP(:,1)+EVFsP qP(:,2) qP(:,3) qP(:,4) qP(:,5)]*1e-6;

phiA  = -cumsum(qA(end:-1:1,:),1);
phiP  = -cumsum(qP(end:-1:1,:),1);

minqA = floor(min(min(phiA)))*ones(n,1) - 2;
minqP = floor(min(min(phiP)))*ones(n,1) - 2;

PsiA = [zeros(n,1) minqA phiA zeros(n,1)];
PsiP = [zeros(n,1) minqP phiP zeros(n,1)];

PsiA = PsiA(end:-1:1,:);
PsiA(end,:) = 0;

PsiP = PsiP(end:-1:1,:);
PsiP(end,:) = 0;

dl = .1*dr;
z=zm*1e-3;

latA = [-Bsh -Bs -fdh-dl -fdh+dl -fdl  fde    fdl   fdh  fBr  Bn]/dr;
latP = [-Bsh -Bs -fdh-dl -fdh+dl -fdl  fde    fdl   fdh  fBr    ]/dr;



xfillAntarc = [-Bsh -pi/2 -pi/2 -Bsh]/dr;
xfillBering = [fBr Bn Bn fBr]/dr;
yfill       = [n     n    0    0]*1e2;


fig=figure('units','normalized','outerposition',[0.1 0.05 .45 .85]);
nf = 3;
nc = 1;
lm  = 0.08;     % left margin
rm  = 0.99;     % right margin
bm  = 0.08;     % bottom margin
tm  = 0.97;     % top margin
hsp = 0.014;    % horizontal space
vsp = 0.03;     % vertical space

dx = rm - lm;
dy = tm - bm;

hd = ( dx-(nc-1)*hsp )/nc; 
vd = ( dy-(nf-1)*vsp )/nf;

k=0;
for i=nf:-1:1
    for j=1:nc
        k=k+1;
        left   = lm + (j-1)*(hd+hsp);
        bottom = bm + (i-1)*(vd+vsp);
        axStream(k,:)  = axes(fig,'Position',[left bottom hd vd],'Box','on');
    end
end

coln = '#00aaff';
colp = '#ff4b4e';


vcp = 1:1:23;
vcn = -21:1:-1;

axes(axStream(1))
[C1,h1]=contour(latA,z,PsiA,vcn,'Color',coln,'Linewidth',1.2); hold on
clabel(C1,h1,vcn(1:2:end));
[C1,h1]=contour(latA,z,PsiA,vcp,'Color',colp,'Linewidth',1.2); hold on
clabel(C1,h1,vcp(1:2:end));
fl=fill(xfillAntarc,yfill,[.5 .5 .5]); set(fl,'edgecolor','none')
axis([-90 90 0 n*1e2*1e-3])
ylabel('Depth (km)','FontSize',15)
axis ij
set(axStream(1),'XTicklabels',{' '},'YTick',0:.5:5.5,'fontsize',12)
axes(axStream(2))
[C1,h1]=contour(latP,z,PsiP,vcn,'Color',coln,'Linewidth',1.2); hold on
clabel(C1,h1,vcn(1:2:end));
[C1,h1]=contour(latP,z,PsiP,vcp,'Color',colp,'Linewidth',1.2); hold on
clabel(C1,h1,vcp(1:2:end));
fl=fill(xfillAntarc,yfill,[.5 .5 .5]); set(fl,'edgecolor','none')
fl=fill(xfillBering,yfill,[.5 .5 .5]); set(fl,'edgecolor','none')
axis([-90 90 0 n*1e2*1e-3])
ylabel('Depth (km)','FontSize',15)
axis ij
set(axStream(2),'XTicklabels',{' '},'YTick',0:.5:5.5,'XTick',-80:20:80,'fontsize',12)

vcp = 1:1:30;
vcn = -25:1:-1;
PsiPG       = PsiP;
PsiPG(:,10) = 0;
latG = latA;
PsiG = PsiA + PsiPG;
axes(axStream(3))
[C1,h1]=contour(latG,z,PsiG,vcn,'Color',coln,'Linewidth',1.2); hold on
clabel(C1,h1,vcn(1:2:end));
[C1,h1]=contour(latG,z,PsiG,vcp,'Color',colp,'Linewidth',1.2); hold on
clabel(C1,h1,vcp(1:2:end));
fl=fill(xfillAntarc,yfill,[.5 .5 .5]); set(fl,'edgecolor','none')
axis([-90 90 0 n*1e2*1e-3])
xlabel('Latitude (degrees)')
ylabel('Depth (km)','FontSize',15)
axis ij
set(axStream(3),'YTick',0:.5:5.5,'fontsize',12)

text(axStream(1) ,.01,.95,'(a)','Units','normalized','FontWeight','bold','FontSize',12);
text(axStream(2) ,.01,.95,'(b)','Units','normalized','FontWeight','bold','FontSize',12);
text(axStream(3) ,.01,.95,'(c)','Units','normalized','FontWeight','bold','FontSize',12);



end