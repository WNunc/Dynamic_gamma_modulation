function plot_hist_shandreal(Qsh,Qreal,bin)
histogram(Qsh,bin)
CI95 = round(prctile(Qsh,95),4);
xline(CI95,'k-','LineWidth',1.5)
xline(Qreal,'r-','LineWidth',1.5)
end