function BorderLine1(obj, dir)
limX = get(obj, 'Xlim');
limY = get(obj, 'Ylim');
minX = limX(1);
maxX = limX(2);
minY = limY(1);
maxY = limY(2);
    if "right" == dir 
        line([maxX,maxX], [minY,maxY],"LineStyle","-","Color", "k","LineWidth",1);
    elseif "left" == dir
        line([minX,minX], [minY,maxY],"LineStyle","-", "Color", "k","LineWidth",1);
    elseif "upper" == dir
        line([minX,maxX], [maxY,maxY],"LineStyle","-", "Color", "k","LineWidth",1);
    end
end