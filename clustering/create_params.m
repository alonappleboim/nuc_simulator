arr = cell(1,100);
intensity = 1:10;
length = 10:10:100;
for i = 0:99
    arr(1,i+1) = {[intensity(1 + (fix(i/10))), length(mod(i,10) + 1)]};
end
        