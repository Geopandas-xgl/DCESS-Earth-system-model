% clear all 

global sy R13pdb mgt dm d n

GAw = load('GAw_Eo.txt');
GAc = load('GAc_Eo.txt');

load Data\senanal_results.mat
load ResThilda_M.mat

warning off
ParVal_M

num_idx = numel(run);
%% 
for i = 1: num_idx
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

[xmin,xmax,tick_min,tick_inter,tick_max] = deal(-50,71,-50,10,71);

[~,end_idx] = min(abs(t_kyr-70.7 )); 
%% 
flux_meth = zeros(num_idx,length(t_kyr));
flux_mett = zeros(num_idx,length(t_kyr));
flux_lip = zeros(num_idx,length(t_kyr));
meth_total_serial = zeros(num_idx,length(t_kyr));
mett_total_serial = zeros(num_idx,length(t_kyr));
lip_total_serial = zeros(num_idx,length(t_kyr));

for i = fliplr(1:num_idx)
    flux_lip(i,:) = 2*run(i).state.scar(4,:)/mgt*sy;
    flux_meth(i,:) = 2*run(i).state.scar(1,:)/mgt*sy;
    flux_mett(i,:) = 2*run(i).state.scar(2,:)/mgt*sy;

    meth_total_serial(i,:) = cumtrapz(run(i).state.st,flux_meth(i,:));
    mett_total_serial(i,:) = cumtrapz(run(i).state.st,flux_mett(i,:));
    lip_total_serial(i,:) = cumtrapz(run(i).state.st,flux_lip(i,:));
end

range.meth = [run(1).pars.meth.t_start,run(1).pars.meth.t_end]/1000 - run(1).pars.time_onset/1000;
range.mett = [run(1).pars.mett.t_start,run(1).pars.mett.t_end]/1000 - run(1).pars.time_onset/1000;
range.lip = [run(1).pars.lip.t_start,run(1).pars.lip.t_end]/1000 - run(1).pars.time_onset/1000;

[~,idx.meth(1)] = min(abs(t_kyr-range.meth(1)));
[~,idx.meth(2)] = min(abs(t_kyr-range.meth(2)));
[~,idx.mett(1)] = min(abs(t_kyr-range.mett(1)));
[~,idx.mett(2)] = min(abs(t_kyr-range.mett(2)));
[~,idx.lip(1)] = min(abs(t_kyr-range.lip(1)));
[~,idx.lip(2)] = min(abs(t_kyr-range.lip(2)));

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

EdgeColor = ["#6A7A5B","#3B6890","#915A59"];

FaceColor = hex2grb(["#DDE1D9","#CDD5DC","#E0D0CF"]);

lnstl = [":","--","-"];


% [mid_idx,~] = middleValue(parsSens.k_A);
mid_idx = 2; % 

figure("Position",[0,0,530,720])
[ha, pos]= tight_subplot(6,1,[0,0],[0.075,0.025],[0.15,0.15]);


axes(ha(1))

% fill([t_kyr(1:idx.meth(2)),fliplr(t_kyr(1:idx.meth(2)))],[flux_lip(1,1:idx.meth(2)),fliplr(flux_lip(num_idx,1:idx.meth(2)))],FaceColor(1,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
% line(t_kyr(1:idx.meth(2)),flux_lip(mid_idx,1:idx.meth(2)),"Color",EdgeColor(1), "LineWidth",1.4);
% 
% fill([t_kyr(idx.meth(1):idx.mett(2)),fliplr(t_kyr(idx.meth(1):idx.mett(2)))],[flux_meth(1,idx.meth(1):idx.mett(2)),fliplr(flux_meth(num_idx,idx.meth(1):idx.mett(2)))],FaceColor(2,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); 
% line(t_kyr(idx.meth(1):idx.mett(2)),flux_meth(mid_idx,idx.meth(1):idx.mett(2)),"Color",EdgeColor(2), "LineWidth",1.4);
% 
% fill([t_kyr(idx.mett(1):end_idx),fliplr(t_kyr(idx.mett(1):end_idx))],[flux_mett(1,idx.mett(1):end_idx),fliplr(flux_mett(num_idx,idx.mett(1):end_idx))],FaceColor(3,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); 
% line(t_kyr(idx.mett(1):end_idx),flux_mett(mid_idx,idx.mett(1):end_idx),"Color",EdgeColor(3), "LineWidth",1.4);

fill([t_kyr,fliplr(t_kyr)],[flux_lip(1,:),fliplr(flux_lip(num_idx,:))],FaceColor(1,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
line(t_kyr,flux_lip(mid_idx,:),"Color",EdgeColor(1), "LineWidth",1.4);

fill([t_kyr,fliplr(t_kyr)],[flux_meth(1,:),fliplr(flux_meth(num_idx,:))],FaceColor(2,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
line(t_kyr,flux_meth(mid_idx,:),"Color",EdgeColor(2), "LineWidth",1.4);

fill([t_kyr,fliplr(t_kyr)],[flux_mett(1,:),fliplr(flux_mett(num_idx,:))],FaceColor(3,:),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
line(t_kyr,flux_mett(mid_idx,:),"Color",EdgeColor(3), "LineWidth",1.4);

box off
xlim([xmin,xmax]);
ylim([-0.02,0.3]);
yticks(0:0.05:0.3);
ticklabels("y",yticks,2,"F");
ylabel(["\fontname{SimSun}碳输入速率","\fontname{Times}(GtC yr^{-1})"],"Interpreter","tex");
yrange = get(gca, 'Ylim');

rectangle_norm([0.58,0.85,0.08,0.08],FaceColor(1,:));
line_norm([0.58,0.58+0.08],[0.85+0.08/2,0.85+0.08/2],EdgeColor(1),1.2);
text_norm(0.58+0.1,0.85+0.03,"\fontname{Times}CO2\_\fontname{SimSun}幔源 \fontname{Times}(-5 ‰)",9);
rectangle_norm([0.58,0.7,0.08,0.08],FaceColor(2,:));
line_norm([0.58,0.58+0.08],[0.7+0.08/2,0.7+0.08/2],EdgeColor(2),1.2);
text_norm(0.58+0.1,0.7+0.03,"\fontname{Times}CH4\_\fontname{SimSun}水合物 \fontname{Times}(-60 ‰)",9);
rectangle_norm([0.58,0.55,0.08,0.08],FaceColor(3,:));
line_norm([0.58,0.58+0.08],[0.55+0.08/2,0.55+0.08/2],EdgeColor(3),1.2);
text_norm(0.58+0.1,0.55+0.03,"\fontname{Times}CH4\_\fontname{SimSun}热变质 \fontname{Times}(-45 ‰)",9);

ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none');
uistack(onset, 'bottom');
% SST = line([-5,-5],yrange,"Color","#AF4184","LineStyle","--","LineWidth",1.4);
% uistack(SST, 'bottom');
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
for i = 1: num_idx
del_13C = (run(i).state.sLL(5,:,:)./run(i).state.sLL(4,:,:)./R13pdb-1)*1e3; 
del_13C = squeeze(del_13C(1,:,:));
del_13C_slw(i,:) = mean(del_13C(1:2,:),1);
end
fill([t_kyr(1:end_idx),fliplr(t_kyr(1:end_idx))],[del_13C_slw(1,1:end_idx),fliplr(del_13C_slw(num_idx,1:end_idx))],hex2grb("#E6E6E6"),"LineWidth",1,"EdgeColor","none","FaceAlpha",0.9); hold on
line(t_kyr(1:end_idx),del_13C_slw(mid_idx,1:end_idx),"Color","#333333", "LineWidth",1.4);
box off
sel_cri = (data.age2 >= xmin+10) & (data.age2 <= xmax);
dC_sel = data.dC2(sel_cri);
Age_sel = data.age2(sel_cri);
plt1 = plot(Age_sel(~isnan(dC_sel)), dC_sel(~isnan(dC_sel)),"Linestyle","-","LineWidth",1.4,"Color","#59A572"); hold on
plt2 = plot(data.age1,data.dC1,"LineWidth",1.4,"Color","#F6BA92"); 
xlim([xmin,xmax]);
ylim([-6,4]);
yticks(-4:2:4);
ticklabels("y",yticks,2,"F");
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
BorderLine1(gca,"left");
set(gca,"YAxisLocation","right","Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");


axes(ha(3))
sAT = zeros(num_idx,length(t_kyr));
for i = 1:num_idx
sAT(i,:) = squeeze(run(i).state.sAT(4,1,:))*1e6/100;
end
fill([t_kyr(1:end_idx),fliplr(t_kyr(1:end_idx))],[sAT(1,1:end_idx),fliplr(sAT(num_idx,1:end_idx))],hex2grb("#E6E6E6"),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
line(t_kyr(1:end_idx),sAT(mid_idx,1:end_idx),"Color","#333333", "LineWidth",1.4);
box off
xlim([xmin,xmax]);
ylim([5,25])
yticks(10:5:25);
yticklabels(["10","15","20","25"]);
ylabel({"  pCO_2", "(× 10^{2} ppm)"})
yrange = get(gca, 'Ylim');
ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none');
uistack(onset, 'bottom'); 
% SST = line([-5,-5],yrange,"Color","#AF4184","LineStyle","--","LineWidth",1.4);
% uistack(SST, 'bottom');
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
text_norm(0.03,0.95,label(3,1),13);
CO2 = line([xmin,xmax],[19.80,19.80],"Color","#E9536A","LineStyle","--","LineWidth",1.4);
uistack(CO2, 'bottom');
text(-35,22,{"pCO_2 in PETM"},"FontName","Times","FontSize",10,"Color","#E9536A");
BorderLine1(gca,"right");
set(gca,"YAxisLocation","left","Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");

axes(ha(4))
yyaxis right;
Tocean = zeros(num_idx,length(t_kyr));
for i = 1:num_idx
Tocean(i,:) = 0.84*squeeze(run(i).state.sLL(1,1,:))+0.16*squeeze(run(i).state.sHL(1,1,:)); % % Surface temperature (<200 m)
end
fill([t_kyr(1:end_idx),fliplr(t_kyr(1:end_idx))],[Tocean(1,1:end_idx),fliplr(Tocean(num_idx,1:end_idx))],hex2grb("#E6E6E6"),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
line(t_kyr(1:end_idx),Tocean(mid_idx,1:end_idx),"Color","#333333", "LineWidth",1.4);
box off
xlim([xmin,xmax]);
ylim([26,34]);
yticks(26:2:34);
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
ylim([-2,6]);
yticks(-2:2:6);
yticklabels(yticks);
ylabel({"\fontname{Times}SST\fontname{SimSun}异常 \fontname{Times}(℃)"},"Interpreter","tex");
BorderLine1(gca,"right");
ax = gca;
ax.YColor = "#E9536A";
% text(-20,3,"Initial warming","Color","#E4DB84","FontName","Times","FontSize",10);
set(gca,"Xcolor","none","Xtick",[],"linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");

axes(ha(5))
mva8 = zeros(num_idx,length(t_kyr));
for i = 1:num_idx
mva8(i,:) = Ocean_glob_mean(8,run(i).state.sLL,run(i).state.sHL)*1000;
end
fill([t_kyr(1:end_idx),fliplr(t_kyr(1:end_idx))],[mva8(1,1:end_idx),fliplr(mva8(num_idx,1:end_idx))],hex2grb("#E6E6E6"),"LineWidth",1,"EdgeColor","none","FaceAlpha",1); hold on
line(t_kyr(1:end_idx),mva8(mid_idx,1:end_idx),"Color","#333333", "LineWidth",1.4);
box off
xlim([xmin,xmax]);
ylim([60,120]);
yticks(60:10:100);
ticklabels("y",yticks,2,"F");
% yticklabels(["60","80","100"]);
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
text(-22,65,"\fontname{SimSun}启动阶段","Color","#FECDD8","FontName","Times","FontSize",10,"Interpreter","tex");
text(7,120,{"   \fontname{SimSun}底栖","\fontname{SimSun}有孔虫灭绝"},"Color","#E88471","FontName","Times","FontSize",10,"Interpreter","tex");
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
for i = 1:num_idx
fanox(i,:) = run(i).state.fanox_global;
end

fill([t_kyr(1:end_idx),fliplr(t_kyr(1:end_idx))],[fanox(1,1:end_idx),fliplr(fanox(num_idx,1:end_idx))],hex2grb("#E6E6E6"),"LineWidth",1,"EdgeColor","none","FaceAlpha",0.8); hold on
line(t_kyr(1:end_idx),fanox(mid_idx,1:end_idx),"Color","#333333", "LineWidth",1.4);

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
text_norm(0.6+0.1,0.7+0.03,"\fontname{Times}U\fontname{SimSun}反演模型",10);
rectangle_norm([0.6,0.55,0.08,0.08],"#E6E6E6");
line_norm([0.6,0.6+0.08],[0.55+0.08/2,0.55+0.08/2],"#333333",1.2);
text_norm(0.6+0.1,0.55+0.03,"\fontname{Times}DCESS\fontname{SimSun}模型",10);

ext = rectangle("Position",[ext_range(1),yrange(1),diff(ext_range),diff(yrange)],'FaceColor', "#E88471", 'EdgeColor', 'none');
uistack(ext,"bottom");
onset = rectangle("Position",[PETM_onset(1),yrange(1),diff(PETM_onset),diff(yrange)],'FaceColor', "#FECDD8", 'EdgeColor', 'none'); hold on
uistack(onset, 'bottom'); 
% POE = rectangle("Position",[POE_range(1),yrange(1),diff(POE_range),diff(yrange)],'FaceColor', "#E4DB84", 'EdgeColor', 'none');
% uistack(POE, 'bottom');
ylabel({"\fontname{SimSun}缺氧海底", "\fontname{Times}(f_{anox}) (%)"},"Interpreter","tex");
text_norm(0.93,0.95,label(6,1),13);
xtickangle(0);
BorderLine1(gca,"left");
xlabel("\fontname{SimSun}相对\fontname{Times}CIE\fontname{SimSun}启动的年龄 \fontname{Times}(kyr)","Interpreter","tex");
set(gca,"YAxisLocation","right","linewidth", 1,"FontSize",10,"FontName", "Times","TickLength",[0.015,0.03],"Layer","top");
%% 
print(gcf,"Figure\DCESS_PETM_1_Chinese","-dpng","-r600");
%% 
% Extract the dissolved oxygen concentration 
ZLL = permute(run(2).state.sLL(8,:,1:end_idx),[2,3,1])*1000; % dissolved oxygen concentration in LL mmol/m3
ZHL = permute(run(2).state.sHL(8,:,1:end_idx),[2,3,1])*1000; % dissolved oxygen concentration in HL mmol/m3

zcent= [dm/2 dm+(d:d:(n-1)*d)-d/2]/1e3;   % Vertical center of boxes
[X, Y] = meshgrid(t_kyr(1:end_idx), zcent);
% xmin = 0;
% xmax = 94;
xtick_inter = 10; 
% onset = -50;

%--- Define unified 10-unit color levels ---
Min_LL = floor(min([ZLL(:)])/10)*10;
Max_LL = ceil(max([ZLL(:)])/10)*10;

Min_HL = floor(min(ZHL(:))/10)*10;
Max_HL = floor(max(ZHL(:))/10)*10;

levels_LL = Min_LL:10:Max_LL;
clim_LL = [levels_LL(1), levels_LL(end)];

levels_HL = Min_HL:10:Max_HL;
clim_HL = [levels_HL(1),levels_HL(end)];

%--- Load colormap (ocean-oxygen) ---
try
    mapOxy = cmocean('oxy', 256);  % Preferred ocean oxygen colormap
catch
    mapOxy = parula(256);          % Fallback
end
%% 
figure("Position",[0,0,840,280])
% [ha,~] = tight_subplot(1,2,[0 0.17],[0.03 0.08],[0.05 0.11]);
[ha,~] = tight_subplot(1,2,[0 0.17],[0.16 0.08],[0.05 0.11]);
axes(ha(1))
contourf(X, Y, ZLL, levels_LL, 'LineColor','none')
colormap(ha(1), mapOxy);
clim(clim_LL);
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
xticklabels(string(tick_min:tick_inter:tick_max));
xtickangle(0);
ylim([0.1,5]);
yticks(0:0.5:5);
ticklabels("y",yticks,2,"F");
ylabel("\fontname{SimSun}深度 \fontname{Times}(km)","Interpreter","tex");
title("\fontname{SimSun}中低纬度区域","Interpreter","tex")
box off
BorderLine(gca,"LineWidth",1,"Side","yre");
set(gca, "YDir","reverse",'TickDir', 'out', "LineWidth",1, "FontSize",10, "FontName","Times", ...
    "TickLength",[0.02,0.03], "Layer","top")
%--- Colorbar outside right of high-latitude subplot ---
cb = colorbar('Position',[0.41 0.1 0.02 0.82]);  % Adjust as needed
cb.Label.String = 'O_2 (mmol/m^3)';
cb.Ticks = clim_LL(1):20:clim_LL(2);
cb.TickLabels = ["20","","60","","100","","140","","180","","220","","260","","300"];
cb.TickDirection = 'out';
cb.FontSize = 9;
xlabel("\fontname{SimSun}相对\fontname{Times}CIE\fontname{SimSun}启动的年龄 \fontname{Times}(kyr)","Interpreter","tex");

%--- High latitude panel ---
axes(ha(2))
contourf(X, Y, ZHL, levels_HL, 'LineColor','none')
colormap(ha(2), mapOxy)
clim(clim_HL)
xlim([xmin,xmax]);
xticks(tick_min:tick_inter:tick_max);
xticklabels(string(tick_min:tick_inter:tick_max));
xtickangle(0);
ylim([0.1,5]);
yticks(0:0.5:5);
ticklabels("y",yticks,2,"F");
ylabel("\fontname{SimSun}深度 \fontname{Times}(km)","Interpreter","tex");
title("\fontname{SimSun}高纬度区域","Interpreter","tex");
box off
BorderLine(gca,"LineWidth",1,"Side","yre");
set(gca, "YDir","reverse", 'TickDir', 'out', "LineWidth",1, "FontSize",10, "FontName","Times", ...
    "TickLength",[0.02,0.03], "Layer","top")

%--- Colorbar outside right of high-latitude subplot ---
cb = colorbar('Position',[0.905 0.1 0.02 0.82]);  % Adjust as needed
cb.Label.String = 'O_2 (mmol/m^3)';
cb.Ticks = clim_HL(1):20:clim_HL(2);
cb.TickLabels = ["80","","120","","160","","200","","240","","280",""];
cb.TickDirection = 'out';
cb.FontSize = 9;
xlabel("\fontname{SimSun}相对\fontname{Times}CIE\fontname{SimSun}启动的年龄 \fontname{Times}(kyr)","Interpreter","tex");

print(gcf,"Figure\Oxygen concentration_Contour_Chinese","-dpng","-r600");

% exportgraphics(gcf, "Figure\Dissolve O2.pdf", 'ContentType', 'vector');
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
%% 
meth_max = max(meth_total_serial,[],2)
format longG
disp(meth_max)
% max(lip_total_serial,1)


%% 
function BorderLine(ax, opts)
    arguments
        ax matlab.graphics.axis.Axes
        opts.LineWidth (1,1) double = 1
        opts.Side (1,1) string {mustBeMember(opts.Side,["","xre","yre"])} = ""
    end

    % 获取极限
    xL = ax.XLim; yL = ax.YLim;
    hold(ax,'on');

    % 默认画右侧和顶部边线
    if opts.Side == ""
        line(ax, [xL(2) xL(2)], yL, 'Color','k','LineWidth',opts.LineWidth);
        line(ax, xL, [yL(2) yL(2)], 'Color','k','LineWidth',opts.LineWidth);
    end

    % left side
    if opts.Side == "xre"
        line(ax, [xL(1) xL(1)], yL, 'Color','k','LineWidth',opts.LineWidth);
        line(ax, xL, [yL(2) yL(2)], 'Color','k','LineWidth',opts.LineWidth);
    end

    % top side
    if opts.Side == "yre"
        line(ax, xL, [yL(1) yL(1)], 'Color','k','LineWidth',opts.LineWidth);
        line(ax, [xL(2) xL(2)], yL, 'Color','k','LineWidth',opts.LineWidth);
    end

    hold(ax,'off');
end

