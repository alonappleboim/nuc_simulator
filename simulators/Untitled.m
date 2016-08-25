FRS2_seq=sequences_structure(1001,:) ;
FRS2_wt = wt_3h(1001,:);
poly_rates=[0,1,2,3,4,5]; 
nuc_sum=run_simulation_from_genome(FRS2_seq,'report',0,'poly_rate',poly_rates(2)); 
