function Oce_results(name)

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

sArc  = ncread(name,'sarc' ); Arc  = squeeze(sArc(:,:,end));
sSOc  = ncread(name,'ssoc' ); SOc  = squeeze(sSOc(:,:,end));
sSlSO = ncread(name,'sslso'); SlSO = squeeze(sSlSO(end,:));


[~,~,~,~,Glob]=OceanGlobalMean(Arc,HLNA,MLNA,LLNA,LLSA,MLSA,...
                                                               HLNP,MLNP,LLNP,LLSP,MLSP,SOc,SlSO);

disp('Global ocean mean values')
disp(['PO4 (mmol m-3): ' num2str(Glob(3)*1e3)])
disp(['O2  (mmol m-3): ' num2str(Glob(9)*1e3)])
disp(['DIC (mol m-3): ' num2str(Glob(5))])
disp(['ALK (mol m-3): ' num2str(Glob(8))])
end

function [mvalAr,mvalAt,mvalPa,mvalSO,mvalG]=OceanGlobalMean(Arc,HLNA,MLNA,LLNA,LLSA,MLSA,...
                                                                    HLNP,MLNP,LLNP,LLSP,MLSP,SOc,SlSO)

for i=[3 5 8 9]
    [mvalAr(i),mvalAt(i),mvalPa(i),mvalSO(i),mvalG(i)]=OceanValues(i,Arc,HLNA,MLNA,LLNA,LLSA,MLSA,...
                                                                               HLNP,MLNP,LLNP,LLSP,MLSP,SOc,SlSO);
end

end

function [mvalAr,mvalAt,mvalPa,mvalSO,mvalG]=OceanValues(l,Arc,HLNA,MLNA,LLNA,LLSA,MLSA,...
                                                               HLNP,MLNP,LLNP,LLSP,MLSP,SOc,SlSO)


% ===================================================================
% Ocean volumes
% global GAA GAP GASO GAAr

GAA  = load('../Input_Data/GAA.dat' );
GAP  = load('../Input_Data/GAP.dat' );
GAAr = load('../Input_Data/GAAr.dat');
GASO = load('../Input_Data/GASO.dat');


global AoAr AonA AomnA AoenA AoesA AomsA shl dm
global AonP AomnP AoenP AoesP AomsP AoSO Aosh
volAr   = dm*AoAr*sum(GAAr);
volAt   = dm*( AonA *sum(GAA(1,:)) + AomnA*sum(GAA(2,:)) +...
               AoenA*sum(GAA(3,:)) + AoesA*sum(GAA(4,:)) + AomsA*sum(GAA(5,:)) );
volPa   = dm*( AonP *sum(GAP(1,:)) + AomnP*sum(GAP(2,:)) +...
               AoenP*sum(GAP(3,:)) + AoesP*sum(GAP(4,:)) + AomsP*sum(GAP(5,:)) );
volSO   = dm*AoSO*sum(GASO);
volSlSO = dm*Aosh*sum(GASO(1:shl));

volG     = volSO + volAt + volPa + volAr + volSlSO;
% ===================================================================
% Ocean values
vArc = Arc(l,:);
vSOc = SOc(l,:);

vHLNA = HLNA(l,:); vMLNA = MLNA(l,:); vLLNA = LLNA(l,:);
vLLSA = LLSA(l,:); vMLSA = MLSA(l,:);

vHLNP = HLNP(l,:); vMLNP = MLNP(l,:); vLLNP = LLNP(l,:);
vLLSP = LLSP(l,:); vMLSP = MLSP(l,:);

Tf = -1.9;
if l==1
    vSlSO  = Tf;
else
    vSlSO  = SlSO(l);
end


% ===================================================================
% Sum ocean values
vtotsumSlSO  = 0;
for i=1:shl
    vsumSlSO     = vSlSO*GASO(i) ;
    vtotsumSlSO  = vtotsumSlSO  + vsumSlSO;
end

vtotsumAr   = sum(vArc.*GAAr);
vtotsumSO   = sum(vSOc.*GASO);

vtotsumHLNA = sum(vHLNA.*GAA(1,:));
vtotsumMLNA = sum(vMLNA.*GAA(2,:));
vtotsumLLNA = sum(vLLNA.*GAA(3,:));
vtotsumLLSA = sum(vLLSA.*GAA(4,:));
vtotsumMLSA = sum(vMLSA.*GAA(5,:));

vtotsumHLNP = sum(vHLNP.*GAP(1,:));
vtotsumMLNP = sum(vMLNP.*GAP(2,:));
vtotsumLLNP = sum(vLLNP.*GAP(3,:));
vtotsumLLSP = sum(vLLSP.*GAP(4,:));
vtotsumMLSP = sum(vMLSP.*GAP(5,:));

% ===================================================================
% Total ocean values per basin
% Arctic
TotovAr = dm*AoAr*vtotsumAr;
% Southern Ocean
TotovSO = dm*AoSO*vtotsumSO ;
% Atlantic
TotovAt = dm*( AonA *vtotsumHLNA + AomnA*vtotsumMLNA + AoenA*vtotsumLLNA + ...
               AoesA*vtotsumLLSA + AomsA*vtotsumMLSA );
% Pacific
TotovPa = dm*( AonP *vtotsumHLNP + AomnP*vtotsumMLNP + AoenP*vtotsumLLNP + ...
               AoesP*vtotsumLLSP + AomsP*vtotsumMLSP );
% Shelf Southern Ocean
TotovSlSO  = dm*Aosh *vtotsumSlSO ;

TotovG = TotovAr + TotovSO + TotovAt + TotovPa + TotovSlSO;

% ===================================================================
% Mean values per basin
mvalAr  = TotovAr/volAr;
mvalSO  = (TotovSO  + TotovSlSO )/(volSO  + volSlSO );
mvalAt  = TotovAt/volAt;
mvalPa  = TotovPa/volPa;

mvalG  = TotovG/volG;

end