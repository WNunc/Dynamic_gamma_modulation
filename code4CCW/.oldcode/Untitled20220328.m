QQ1 = [];QQ2 = [];
for nl = 1:5
    QQ1 = [QQ1; Q{1,nl}(:,1)];
    QQ2 = [QQ2; Q{1,nl}(:,2)];
end
minq1 = min(QQ1);
minq2 = min(QQ2);
maxq1 = max(QQ1);
maxq2 = max(QQ2);

for nl = 1:5
    Q{1,nl}(:,1) = 2*((Q{1,nl}(:,1)-minq1)/(maxq1-minq1))-1;
    Q{1,nl}(:,2) = 2*((Q{1,nl}(:,2)-minq1)/(maxq1-minq1))-1;
end