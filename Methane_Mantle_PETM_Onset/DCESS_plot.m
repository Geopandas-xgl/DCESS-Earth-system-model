% clear all 

global sy R13pdb mgt dm aLL aHL

GAw = load('GAw_Eo.txt');
GAc = load('GAc_Eo.txt');

load Data\senanal_results.mat
load ResThilda_M.mat

warning off
ParVal_M

num_idx = numel(run);
%% 
for i = [2,4,6]
    [fanox_global, omz_global] = DCESS_fanox(run(i));
    run(i).state.fanox_global = fanox_global;
    run(i).state.omz_global = omz_global;
end
%% 
data1 = readtable("Data\Xiong-2023-PETM_1.xlsx","Sheet","C and O isotope");
data2 = readtable("Data\Site 690.xlsx","Sheet","site 690_Westerhold");
data3 = readtable("Data\Xiong-2025-PETM.xlsx","Sheet","Anoxic area");
data4 = readtable("Data\PETM_0.9_loess.xlsx");
color = readtable("Data\Xiong-2025-PETM.xlsx","Sheet","Color");
recovered_trends = readtable("Data\recovered_trends_Tingri.csv");


data.dC1 = data1.dC(17:end-3);
data.age1 = data1.age_westerhold(17:end-3);   
data.dC2 = data2.d13C;
data.age2 = data2.Age_westerhold;

data.Ba = data4;
data.Ba.series = data.Ba.series/1000 - 22.3;

data.age3 = data3.time/1000;
data.dU3 = data3(:,3:102);
data.fanox3 = data3(:,203:end);

numericCols3 = varfun(@isnumeric,data.fanox3,"OutputFormat","uniform");
data.fanox3{:,numericCols3} = 100*10.^data.fanox3{:,numericCols3};

data.dU_age3 = [[data.age3,data.dU3.d238U_0_16];[flip(data.age3),flip(data.dU3.d238U_0_84)]];
data.fanox_age3 = [[data.age3,data.fanox3.fanox_0_16];[flip(data.age3),flip(data.fanox3.fanox_0_84)]];
data.num_idx3 = 1:size(data.dU_age3,1);

t_kyr = run(1).state.st/1e3 - run(1).pars.time_onset/1e3;

age_initial = -52.77; 
time_seq = table2array(recovered_trends(:,2))/1000 - 22.3;

d238Usw_quant_out = table2array(recovered_trends(:,3:102));
Nsw_quant_out = table2array(recovered_trends(:,103:202));
fanox_quant_out = table2array(recovered_trends(:,203:302));

% [xmin,xmax,tick_min,tick_inter,tick_max] = deal(-60,200,-50,50,200);
[xmin,xmax,tick_min,tick_inter,tick_max] = deal(-50,72,-50,10,72);

[~,end_idx] = min(abs(t_kyr-190)); 
%% 
flux_meth = zeros(num_idx,length(t_kyr));
flux_mett = zeros(num_idx,length(t_kyr));
flux_mix = zeros(num_idx,length(t_kyr));
flux_lip = zeros(num_idx,length(t_kyr));
meth_total_serial = zeros(num_idx,length(t_kyr));
mett_total_serial = zeros(num_idx,length(t_kyr));
lip_total_serial = zeros(num_idx,length(t_kyr));
mix_total_serial = zeros(num_idx,length(t_kyr));

for i = 1:num_idx
    flux_lip(i,:) = 2*run(i).state.scar(4,:)/mgt*sy;
    flux_mix(i,:) = 2*run(i).state.scar(3,:)/mgt*sy;
    flux_meth(i,:) = 2*run(i).state.scar(1,:)/mgt*sy;
    flux_mett(i,:) = 2*run(i).state.scar(2,:)/mgt*sy;
    
    meth_total_serial(i,:) = cumtrapz(run(i).state.st,flux_meth(i,:));
    mett_total_serial(i,:) = cumtrapz(run(i).state.st,flux_mett(i,:));
    lip_total_serial(i,:) = cumtrapz(run(i).state.st,flux_lip(i,:));
    mix_total_serial(i,:) = cumtrapz(run(i).state.st,flux_mix(i,:));
end

range.meth = [run(1).pars.meth.t_start,run(1).pars.meth.t_end]/1000 - run(1).pars.time_onset/1000;
range.mett = [run(1).pars.mett.t_start,run(1).pars.mett.t_end]/1000 - run(1).pars.time_onset/1000;
range.lip = [run(1).pars.lip.t_start,run(1).pars.lip.t_end]/1000 - run(1).pars.time_onset/1000;
range.mix = [run(1).pars.mix.t_start,run(1).pars.mix.t_end]/1000 - run(1).pars.time_onset/1000;

[~,idx.meth(1)] = min(abs(t_kyr-range.meth(1)));
[~,idx.meth(2)] = min(abs(t_kyr-range.meth(2)));
[~,idx.mett(1)] = min(abs(t_kyr-range.mett(1)));
[~,idx.mett(2)] = min(abs(t_kyr-range.mett(2)));
[~,idx.lip(1)] = min(abs(t_kyr-range.lip(1)));
[~,idx.lip(2)] = min(abs(t_kyr-range.lip(2)));
[~,idx.mix(1)] = min(abs(t_kyr-range.mix(1)));
[~,idx.mix(2)] = min(abs(t_kyr-range.mix(2)));

% idx.meth = find(t_kyr(range.meth))

%% 
label = ["A","B","C","D","E","F"];
row = 1:6;
col = 1;
label = reshape(label,length(row),length(col));

PETM_onset = [-4, 5];
SST_onset = [-5,-4];
POE_range = [-37.5,-35.5];
ext_range = [1,4]; 


EdgeColor = ["#6A7A5B","#3B6890","#915A59","#F6BA92"];

FaceColor = hex2grb(["#DDE1D9","#CDD5DC","#E0D0CF","#FFDBB7"]);

% lnstl = [" "," ","-","--","-.",":"];
lnstl(2) = ":";
lnstl(4) = "-";
lnstl(6) = "--";


mid_idx = 4;

figure("Position",[0,0,530,720])
[ha, pos]= tight_subplot(6,1,[0,0],[0.075,0.025],[0.15,0.15]);

axes(ha(1))
line(t_kyr,flux_lip(mid_idx,:),"Color",EdgeColor(1), "LineWidth",1.4); hold on

% fill([t_kyr,fliplr(t_kyr)],[flux_mix(1,:),fliplr(flux_mix(num_idx,:))],FaceColor(1,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
% line(t_kyr,flux_mix(mid_idx,:),"Color",EdgeColor(1), "LineWidth",1.4);
% 
% fill([t_kyr,fliplr(t_kyr)],[flux_meth(1,:),fliplr(flux_meth(num_idx,:))],FaceColor(2,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
% line(t_kyr,flux_meth(mid_idx,:),"Color",EdgeColor(2), "LineWidth",1.4);

line(t_kyr,flux_mett(mid_idx,:),"Color",EdgeColor(3), "LineWidth",1.4);


for i = [2,4,6]
    line(t_kyr(idx.mix(1):end_idx),flux_mix(i,idx.mix(1):end_idx),"Color",EdgeColor(1), "LineWidth",1.4,"LineStyle",lnstl(i));
    line(t_kyr(idx.meth(1):idx.mett(2)),flux_meth(i,idx.meth(1):idx.mett(2)),"Color",EdgeColor(2), "LineWidth",1.4,"LineStyle",lnstl(i));
end
box off
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
ylim([-0.02,0.7])
yticks(0:0.1:0.7);
ticklabels('y',yticks,2,"F");
ylabel(["Carbon input rate","   (GtC yr^{-1})"]);
yrange = get(gca, 'Ylim');
% rectangle_norm([0.6,0.85,0.08,0.08],FaceColor(1,:));
line_norm([0.6,0.6+0.08],[0.85+0.08/2,0.85+0.08/2],EdgeColor(1),1.2);
text_norm(0.6+0.1,0.85+0.03,"CO2\_mantle (-5 ‰)",9);
% rectangle_norm([0.6,0.7,0.08,0.08],FaceColor(2,:));
line_norm([0.6,0.6+0.08],[0.7+0.08/2,0.7+0.08/2],EdgeColor(2),1.2);
text_norm(0.6+0.1,0.7+0.03,"CH4\_hydrat (-60 ‰)",9);
% rectangle_norm([0.6,0.55,0.08,0.08],FaceColor(3,:));
line_norm([0.6,0.6+0.08],[0.55+0.08/2,0.55+0.08/2],EdgeColor(3),1.2);
text_norm(0.6+0.1,0.55+0.03,"CH4\_therm (-45 ‰)",9);
% rectangle_norm([0.6,0.40,0.08,0.08],FaceColor(4,:));
% line_norm([0.6,0.6+0.08],[0.4+0.08/2,0.4+0.08/2],EdgeColor(1),1.2);
% text_norm(0.6+0.1,0.4+0.03,"CO2\_ (-11 ‰)",9);

ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none');
uistack(onset, 'bottom');
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
text_norm(0.03,0.9,label(1,1),13);
lgd.ItemTokenSize = [18,1];
lgd.BackgroundAlpha = 0;
BorderLine1(gca,"upper");
BorderLine1(gca,"right");
set(gca,"YAxisLocation","left","Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");

axes(ha(2))
del_13C_slw = zeros(num_idx,length(t_kyr));
for i = [2,4,6]
    del_13C = (run(i).state.sLL(5,:,:)./run(i).state.sLL(4,:,:)./R13pdb-1)*1e3; 
    del_13C = squeeze(del_13C(1,:,:));
    del_13C_slw(i,:) = mean(del_13C(1:2,:),1);
    line(t_kyr(1:end_idx),del_13C_slw(i,1:end_idx),"Color","#333333", "LineWidth",1.4,"LineStyle",lnstl(i)); hold on
end
box off
sel_cri = (data.age2 >= xmin+10) & (data.age2 <= xmax);
dC_sel = data.dC2(sel_cri);
Age_sel = data.age2(sel_cri);
% plt1 = plot(Age_sel(~isnan(dC_sel)), dC_sel(~isnan(dC_sel)),"Linestyle","-","LineWidth",1.4,"Color","#59A572"); hold on
% plt2 = plot(data.age1,data.dC1,"LineWidth",1.4,"Color","#F6BA92"); 
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
ylim([-6,4]);
yticks(-4:2:4);
yticklabels(yticks);
% ticklabels("y",yticks,2,"F");
ylabel('\delta^{13}C (‰)');
yrange = get(gca, 'Ylim');
ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none');
uistack(onset, 'bottom'); 
% SST = line([-5,-5],yrange,"Color","#AF4184","LineStyle","--","LineWidth",1.4);
% uistack(SST, 'bottom');
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
text_norm(0.93,0.9,label(2,1),13);
text(7,-3,{"  Benthic","extinction"},"Color","#E88471","FontName","Times","FontSize",10);
BorderLine1(gca,"left");
set(gca,"YAxisLocation","right","Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");


axes(ha(3))
sAT = zeros(num_idx,length(t_kyr));
for i = [2,4,6]
    sAT(i,:) = squeeze(run(i).state.sAT(4,1,:))*1e6/100;
    line(t_kyr(1:end_idx),sAT(i,1:end_idx),"Color","#333333", "LineWidth",1.4,"LineStyle",lnstl(i));
end
% fill([t_kyr(1:end_idx),fliplr(t_kyr(1:end_idx))],[sAT(1,1:end_idx),fliplr(sAT(num_idx,1:end_idx))],hex2grb("#E6E6E6"),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
% line(t_kyr(1:end_idx),sAT(mid_idx,1:end_idx),"Color","#333333", "LineWidth",1.4);
box off
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
ylim([5,27])
yticks(10:5:25);
yticklabels(["10","15","20","25"]);
ylabel({"  pCO_2", "(× 10^{2} ppm)"})
yrange = get(gca, 'Ylim');
xrange = get(gca, "Xlim");
ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none');
uistack(onset, 'bottom'); 
% SST = line([-5,-5],yrange,"Color","#AF4184","LineStyle","--","LineWidth",1.4);
% uistack(SST, 'bottom');
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
CO2 = line(xrange,[19.80,19.80],"Color","#E9536A","LineStyle","--","LineWidth",1.4);
uistack(CO2, 'bottom');
text(-35,22,{"pCO_2 in PETM"},"FontName","Times","FontSize",10,"Color","#E9536A");
text_norm(0.03,0.95,label(3,1),13);
BorderLine1(gca,"right");
set(gca,"YAxisLocation","left","Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");

axes(ha(4))
yyaxis right;
Tocean = zeros(num_idx,length(t_kyr));
for i = [2,4,6]
    Tocean(i,:) = 0.84*squeeze(run(i).state.sLL(1,1,:))+0.16*squeeze(run(i).state.sHL(1,1,:)); % % Surface temperature (<200 m)
    line(t_kyr(1:end_idx),Tocean(i,1:end_idx),"Color","#333333", "LineWidth",1.4,"LineStyle",lnstl(i)); hold on
end 
box off
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
ylim([26,36]);
yticks(26:2:36);
yticklabels(yticks);
ylabel('SST (^oC)');
yrange = get(gca, 'Ylim');
ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none');
uistack(onset, 'bottom'); 
% SST = line([-5,-5],yrange,"Color","#AF4184","LineStyle","--","LineWidth",1.4);
% uistack(SST, 'bottom');
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
text_norm(0.93,0.95,label(4,1),13);
xtickangle(0);
BorderLine1(gca,"left");
ax = gca;
ax.YColor = "#000000";
set(gca,"Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");

yyaxis left;
box off
xlim([xmin,xmax]);
ylim([-2,8]);
yticks(-2:2:8);
yticklabels(yticks);
ylabel({"SST anomaly", "   (℃)"});
BorderLine1(gca,"right");
ax = gca;
ax.YColor = "#E9536A";
set(gca,"Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");

axes(ha(5))
mva8 = zeros(num_idx,length(t_kyr));
for i = [2,4,6]
    mva8(i,:) = Ocean_glob_mean(8,run(i).state.sLL,run(i).state.sHL)*1000;
    line(t_kyr(1:end_idx),mva8(i,1:end_idx),"Color","#333333", "LineWidth",1.4,"LineStyle",lnstl(i));
end
box off
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
ylim([60,120]);
yticks(60:10:100);
ticklabels("y",yticks,2,"F");
ylabel({"O_2", "(mmol m^{-3})"});
yrange = get(gca, 'Ylim');
ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none');
uistack(onset, 'bottom'); 
% SST = line([-5,-5],yrange,"Color","#AF4184","LineStyle","--","LineWidth",1.4);
% uistack(SST, 'bottom');
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
text_norm(0.03,0.95,label(5,1),13);
% text(-33,65,"POE","Color","#E4DB84","FontName","Times","FontSize",10);
text(-15,65,"Onset","Color","#FECDD8","FontName","Times","FontSize",10);
BorderLine1(gca,"right");
set(gca,"YAxisLocation","left","Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");

axes(ha(6))
for i = 1:45
    lower = 100*10.^fanox_quant_out(:, 50 - i + 1);
    upper = 100*10.^fanox_quant_out(:, 50 + i);
    patch([time_seq; flipud(time_seq)], [lower; flipud(upper)], [55 153 181]/255,'FaceAlpha', 0.1,'EdgeColor', 'none'); hold on
end
lower16 = 100*10.^fanox_quant_out(:,16);
upper84 = 100*10.^fanox_quant_out(:,84);
patch([time_seq; flipud(time_seq)],[lower16; flipud(upper84)],[55 153 181]/255,"FaceColor",'none','EdgeColor', [55 153 181]/255,'LineStyle', '--','LineWidth', 1); 
median50 = 100*10.^fanox_quant_out(:,50);
plot(time_seq, median50, 'Color', '#006A87', 'LineWidth', 1.2);
fanox = zeros(num_idx,length(t_kyr));
for i = [2,4,6]
    fanox(i,:) = run(i).state.fanox_global;
    line(t_kyr(1:end_idx),fanox(i,1:end_idx),"Color","#333333", "LineWidth",1.4,"LineStyle",lnstl(i));
end
box off
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
ticklabels("x",xticks,2,"T");
ylim([0,5]);
yticks(0:1:5);
yticklabels(yticks);
yrange = get(gca, 'Ylim');

rectangle_norm([0.6,0.7,0.08,0.08],"#A0CFDC");
line_norm([0.6,0.6+0.08],[0.7+0.08/2,0.7+0.08/2],'#006A87',1.2);
text_norm(0.6+0.1,0.7+0.03,"U inverse model",10);
% rectangle_norm([0.6,0.55,0.08,0.08],"#E6E6E6");
% line_norm([0.6,0.6+0.08],[0.55+0.08/2,0.55+0.08/2],"#333333",1.2);
% text_norm(0.6+0.1,0.55+0.03,"DCESS model",10);
ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none'); hold on
uistack(onset, 'bottom'); 
% SST = line([-5,-5],yrange,"Color","#AF4184","LineStyle","--","LineWidth",1.4);
% uistack(SST, 'bottom');
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
ylabel({"Anoxic seafloor", "(f_{anox}) (%)"});
% ylabel("Seafloor anoxic (%)");
text_norm(0.93,0.95,label(6,1),13);
xtickangle(0);
BorderLine1(gca,"left");
xlabel("Age relative to CIE onset (kyr)","FontSize",14);
set(gca,"YAxisLocation","right","linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");
%% 
print(gcf,"Figure\DCESS_PETM","-dpng","-r600");
%% 
[~,idx1] = min(abs(t_kyr-22.3));
C1 = run(1).state.sAT(4,1,idx1)*1e6;

T1 = Tocean(idx1);

C2 = max(run(1).state.sAT(4,1,:))*1e6;
T2 = max(Tocean);

dC = C2-C1;

dT = T2-T1;

[~, idx2] = min(abs(t_kyr -20));

O1 = (max(mva8)-mva8(idx2))/max(mva8)*100;

max(mix_total_serial,[],2)
