



earlyPhs = [0.0898364	0.034566
0.0242126	0.0239977
0.0041209	0.0875909
0.1420617	-0.0106
0.1355113	0.170936
0.1344464	0.1474726
0.049214	0.111975
0.0239314	0.0017711
0.0818975	0.0670467
0.155323	0.1278995
0.004204	0.023722
0.0364639	0.0144906
0.014209	-0.04224
0.126532	0.056154
0.096508	0.089094
0.170692	0.090001
0.126207	0.075697
0.030065	0.198842
0.090764	0.130691
0.141931	0.160231
0.013487	0.111519
0.073885	0.105102
0.154904	0.123923
0.080731	0.036958];

color = {'#F24444','#F2CA50'};
figure('Position', [1179 316 283 278])
for x = 1:2
[bar,er] = barwitherror(x,earlyPhs(:,x));
bar.FaceColor = color{x};
hold on
end
hold off
xticks([1,2])
xticklabels({'FG-cell','NFG-cell'})
set(gca,'FontSize',14)
ylim([0,0.2])
ylabel('Weight correlation')


latePhs = [-0.01431	0.123628
0.128529	0.152775
0.041232	0.276337
0.20621	0.142508
0.126292	0.176064
0.096265	0.193924
0.082073	0.05876
0.111924	0.238598
0.024391	0.02259
0.273467	0.181634
0.197792	0.245407
0.058854	0.105321
0.011529	0.142435
0.154163	0.166447
0.011982	0.141185
0.077695	0.177005
0.18544	0.289778
0.080548	0.114664
0.086323	0.097718
0.069469	0.167332
0.065332	0.057254
0.196604	0.202904
0.055728	0.066306
-0.00218	0.182087];

color = {'#F24444','#F2CA50'};
figure('Position', [1179 316 283 278])
for x = 1:2
[bar,er] = barwitherror(x,latePhs(:,x));
bar.FaceColor = color{x};
hold on
end
hold off
xticks([1,2])
xticklabels({'FG-cell','NFG-cell'})
set(gca,'FontSize',14)
ylim([0,0.2])
ylabel('Weight correlation')
