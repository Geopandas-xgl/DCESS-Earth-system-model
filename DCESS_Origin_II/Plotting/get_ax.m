function ax=get_ax

fig = figure('units','normalized','outerposition',[0.1 0.05 .45 .85]);
nf  = 3;
nc  = 4;
lm  = 0.05;     % left margin
rm  = 0.98;     % right margin
bm  = 0.06;     % bottom margin
tm  = 0.94;     % top margin
hsp = 0.014;    % horizontal space
vsp = 0.03;    % vertical space

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
        ax(k,:)    = axes(fig,'Position',[left bottom hd vd],'Box','on');
    end
end



end