function mval=Ocean_glob_mean(l,sLL,sHL)
global dm aLL aHL

GAw = load('GAw_Eo.txt');
GAc = load('GAc_Eo.txt');
GAw = GAw/GAw(1);
GAc = GAc/GAc(1);

Totvol=2*dm*(aLL*sum(GAw)+aHL*sum(GAc));

for i=1:length(sLL(1,1,:))
    valLL=squeeze(sLL(l,:,i));
    valHL=squeeze(sHL(l,:,i));

    sumLL=sum(valLL.*GAw);
    sumHL=sum(valHL.*GAc);
    Totocean(i)= 2*dm*(aLL*sumLL+aHL*sumHL);
    mval(i)    = Totocean(i)/Totvol; % global mean

end

end