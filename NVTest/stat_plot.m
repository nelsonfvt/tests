function [] = stat_plot(Mx)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
figure()
xticks([1 2 3])
xticklabels({'S1','S2','S3'})
xlabel('S-Statistics')
ylabel('Normalized values')
hold on
for i=1:10
    plot(Mx(i,:),'--s')
end
hold off
end

