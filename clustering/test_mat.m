b = zeros(1,10);
for i=1:10
	load(['mat' num2str(i) '.mat'],'a');
	b(i) = a;
end

save('matsum.mat','b');