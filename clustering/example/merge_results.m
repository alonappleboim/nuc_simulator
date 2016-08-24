ab_idx = 1:36;
first = true;
for i = 1:67
    fname = sprintf('nuc_fit_batch_%i.mat',i);
    try
        load(fname)
        if first
            fits = struct();
            fnames = fieldnames(fit);
            vals = cell(1,length(fnames)-1);
            first = false;
        end
        for fi = 1:length(fnames)-1 %last field is indices
            f = fnames{fi};
            tmp = fit.(f);
            vals{fi}(fit.nuc_inds,ab_idx) = tmp(:, ab_idx);
        end
    catch
        fprintf('couldnt read file %s\n', fname)
    end
end

fit = struct();
for fi = 1:length(fnames)-1 %last field is indices
    fit.(fnames{fi}) = vals{fi};
end
clear fits fname fnames tmp i fi f first ab_idx vals;