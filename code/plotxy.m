function plotxy(x,y, fnum, ftitle, xtitle, ytitle, opt_grid, opt_hold, ...
    splotx, sploty)


if splotx > 1 || sploty > 1
    subplot(splotx,sploty,fnum)
else
    figure(fnum)
end


plot(x,y)
title(ftitle)
xlabel(xtitle)
ylabel(ytitle)
grid(opt_grid)
hold(opt_hold)

