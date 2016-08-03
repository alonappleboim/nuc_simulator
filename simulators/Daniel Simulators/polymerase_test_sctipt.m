poly_params = 1:10;
distances = zeros(1,10);
for i=1:10
    e = simulate_polymerase(poly_params(i));
    distances(i) = e;
end

bar(poly_params,distances)