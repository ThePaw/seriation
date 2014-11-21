%% Octave

set (0, "defaultaxesfontname", "Arial")
set (0, "defaulttextfontname", "Arial")


load('HRRankMatrix50.dat')
colormap(copper)
imagesc(HRRankMatrix50);
axis ("square");
caxis([0, 3])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H rank matrix, n=50', 'fontsize',18);
print -dpdf HRRankMatrix50.pdf

load('HRRankMatrix100.dat')
colormap(copper)
imagesc(HRRankMatrix100);
axis ("square");
caxis([0, 3])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H rank matrix, n=100', 'fontsize',18);
print -dpdf HRRankMatrix100.pdf

load('HRRankMatrix150.dat')
colormap(copper)
imagesc(HRRankMatrix150);
axis ("square");
caxis([0, 3])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H rank matrix, n=150', 'fontsize',18);
print -dpdf HRRankMatrix150.pdf

load('HRRankMatrix250.dat')
colormap(copper)
imagesc(HRRankMatrix250);
axis ("square");
caxis([0, 3])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H rank matrix, n=250', 'fontsize',18);
print -dpdf HRRankMatrix250.pdf



load('HRRankMatrix500.dat')
set (0, "defaultaxesfontname", "Arial")
set (0, "defaulttextfontname", "Arial")
colormap(copper)
imagesc(HRRankMatrix500);
axis ("square");
caxis([0, 3])
set(gca,'fontsize',12); % sets font of numbers on axes 
colorbar ("EastOutside")
title('H rank matrix, n=500', 'fontsize',18);
print -dpdf HRRankMatrix500.pdf
