tmp=all.sth1_0(1001,500:1500);

cov=ksdensity([1:length(tmp)],[1:length(tmp)],'weights',double(tmp),'width',30);


[vals_c,peaks_c] =findpeaks(cov,'MinPeakHeight',4*10^-4);

