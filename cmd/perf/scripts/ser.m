%% Octave

set (0, "defaultaxesfontname", "Arial")
set (0, "defaulttextfontname", "Arial")


load('Pre150.dat')
colormap(copper)
imagesc(Pre150);
axis ("square");
%caxis([0, 1])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('Pre-Robinson matrix, n=150', 'fontsize',18);
print -dpdf Pre150.pdf

load('H150.dat')
colormap(copper)
imagesc(H150);
axis ("square");
%caxis([0, 1])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('H-seriated Pre-Robinson matrix, n=150', 'fontsize',18);
print -dpdf H150.pdf

load('R150.dat')
colormap(copper)
imagesc(R150);
axis ("square");
%caxis([0, 1])
set(gca,'fontsize',12); % sets font of numbers on axes 
title('Robinson matrix, n=150', 'fontsize',18);
colorbar ("EastOutside")
print -dpdf R150.pdf

