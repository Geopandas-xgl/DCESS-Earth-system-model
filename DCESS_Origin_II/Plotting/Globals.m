function Globals

global dr sy mgt fdh fdl fde Bsh Bs fBr Bn n zm shl

dr      =  pi/180;      % Degrees to radians
sy      =  31536000;	% Seconds per year
mgt     =  8.326e13;        % Moles C per GtC

fde     =  0*dr; 
fdl     =  35*dr;       % Dividing Latitude, both north and south
fdh     =  55*dr;       % Dividing Latitude, both north and south
fBr     =  65*dr;       % Dividing Latitude, both north and south
Bn   	=  90*dr;		% Northern Ocean boundary
Bs      =  69*dr;       % Southern Ocean boundary
Bsh   	=  70*dr;       % Southern Column boundary

n       =  55;		% Total number of vertical layers
d       =  100;		% Layer thickness 
zm       = (0:n-1)*d+d/2; % centered model layer depths

global dm
dm  = 100;          % Mixed layer depth (m)
shl = 5;            % Shelf level /Inflow- outflow (-)

% Ocean areas
global AonA AomnA AoenA AoesA AomsA
global AonP AomnP AoenP AoesP AomsP
global AoAr AoSO Aosh
% Ocean to land fractions 
olfnA  =  0.2419;
olfmnA =  0.2227;
olfenA =  0.1939;
olfesA =  0.1560;
olfmsA =  0.2338;

olfnP  =  0.1397;
olfmnP =  0.2773;
olfenP =  0.5024;
olfesP =  0.6233;
olfmsP =  0.7289;

olfAr  =  0.5914;

r       =  6.371e6;     % Earth Radius
Ah	    =  2*pi*r^2;    % Hemispheric Area

AonA     =  olfnA *Ah*(sin(fBr) - sin(fdh));
AomnA    =  olfmnA*Ah*(sin(fdh) - sin(fdl));
AoenA    =  olfenA*Ah*(sin(fdl) - sin(fde));
AoesA    =  olfesA*Ah*(sin(fdl) - sin(fde));
AomsA    =  olfmsA*Ah*(sin(fdh) - sin(fdl));


AonP     =  olfnP *Ah*(sin(fBr) - sin(fdh));
AomnP    =  olfmnP*Ah*(sin(fdh) - sin(fdl));
AoenP    =  olfenP*Ah*(sin(fdl) - sin(fde));
AoesP    =  olfesP*Ah*(sin(fdl) - sin(fde));
AomsP    =  olfmsP*Ah*(sin(fdh) - sin(fdl));

AoAr =  olfAr *Ah*(sin(Bn) - sin(fBr));
AoSO = Ah*(sin(Bs)  - sin(fdh));
Aosh = Ah*(sin(Bsh) - sin(Bs) );

global R13pdb R14oas
R13pdb  = 0.0112372;                  % 'Pee-Dee Belemnite' 13C/12C standard []
R14oas  = 1.176e-12;                  % 'Oxalic Acid Standard' for the 14C/12C ratio []
end