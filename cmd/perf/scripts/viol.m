%% Octave

set (0, "defaultaxesfontname", "Arial")
set (0, "defaulttextfontname", "Arial")


load('HRViolMatrix150.dat')
colormap(copper)
imagesc(HRViolMatrix150);
axis ("square");
%caxis([0, 20])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H pair-order violations matrix, n=150', 'fontsize',18);
colorbar ("EastOutside")
print -dpdf HRViolMatrix150.pdf

load('HRViolMatrix250.dat')
colormap(copper)
imagesc(HRViolMatrix250);
axis ("square");
%caxis([0, 20])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H pair-order violations matrix, n=250', 'fontsize',18);
colorbar ("EastOutside")
print -dpdf HRViolMatrix250.pdf

load('HRViolMatrix500.dat')
colormap(copper)
imagesc(HRViolMatrix500);
axis ("square");
%caxis([0, 20])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H pair-order violations matrix, n=500', 'fontsize',18);
colorbar ("EastOutside")
print -dpdf HRViolMatrix500.pdf
