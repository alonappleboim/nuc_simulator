function [ buffered_gene ] = create_gene_buffer( gene, genlen )
%create_gene_buffer create a buffered gene from a 2501-bp gene
%   Given the gene to buffer and the gene length we want, the function
%   returns a new gene vector in the desired length, with the TSS in the
%   middle of the gene. 

buffer = genlen - 2501;
right_buffer = fix((buffer-500)/2);
left_buffer = right_buffer + 500;
if (right_buffer + left_buffer < buffer)
    left_buffer = left_buffer + 1;
end

if (right_buffer < 1 || left_buffer < 1)
    buffered_gene = gene;
else
    buffered_gene = [zeros(1,left_buffer), gene, zeros(1, right_buffer)];
end

end

