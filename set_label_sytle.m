function [] = set_label_sytle(handle, labelFontSize, labelFontWeight, fullScreen)
switch nargin
    case 1
        labelFontSize = 14;
        labelFontWeight = 'normal';
        fullScreen = false;
    case 2
        labelFontWeight = 'normal';
        fullScreen = false;
    case 3
        fullScreen = false;
end

if isempty(labelFontWeight)
    labelFontWeight = 'normal';
end

% Set label style
figure(handle)
set(findall(handle,'type','text'),'fontSize', labelFontSize,'fontWeight', labelFontWeight)
set(gca,'FontSize', labelFontSize, 'fontWeight', labelFontWeight)

% Make figure full screen
if fullScreen
    set(handle, 'units','normalized','outerposition',[0 0 1 1]);
end

end