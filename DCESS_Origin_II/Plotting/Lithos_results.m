function Lithos_results

global sy mgt

% BorgPI  = 2.4830e+03;  % global PO4 burial
BcarPI  = 5.6722e+05;  % global carbonate burial
BCorgPI = 2.7555e+05;  % global organic carbon burial

gamma_sil  = 0.85;
d13VolPI   = -0.005;
d13CcarPI  = 2.2534*1e-3;
d13CorgCPI = -22.0357*1e-3;

% Volcan input
Volo  = gamma_sil/(1+gamma_sil)*BcarPI*(d13CcarPI - d13CorgCPI)/(d13VolPI - d13CorgCPI);
Wcarb = (BcarPI  /(1+gamma_sil));
Worg  = (BCorgPI + gamma_sil*BcarPI/(1+gamma_sil)) - Volo;
Wsil  = gamma_sil*Wcarb;


disp(['Vol Outgassing (Gt C/yr): ' num2str(Volo*sy/mgt)])
disp(['Organic Carbon Weathering (Gt C/yr): ' num2str(Worg*sy/mgt)])
disp(['Carbonate Weathering (Gt C/yr): ' num2str(Wcarb*sy/mgt)])
disp(['Silicate Weathering (Gt C/yr): ' num2str(Wsil*sy/mgt)])

end