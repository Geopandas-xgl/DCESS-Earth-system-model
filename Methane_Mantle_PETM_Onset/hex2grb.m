function rgbColors = hex2grb(hexColors)
rgbColors = zeros(length(hexColors), 3);
for i = 1:length(hexColors)
    rgbColors(i, :) = sscanf(hexColors(i).char, '#%2x%2x%2x', [1 3]) / 255;
end
end