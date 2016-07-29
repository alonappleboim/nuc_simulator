function [ dp ] = ode_i_o( t, p )
%ode_i_o is the ode function for the mode that only has input and output of
%nucleosomes. this function uses non linear differential equations for the
%actual numeric process, to avoid exponential states...

    % nucleosome length for model:
    nuc_len = 100;

    % vectors for the input rate (alpha) and output rate (beta) of each bp
    % in the genome (we assume these are constant throughout the process):
    alpha = 0.1 * ones(size(p));
    beta  = 0.1 * ones(size(p));
    
    % don't allow entry of nucleosomes in the edges:
    for i=1:int8(nuc_len/2)
        alpha(i) = 0;
        alpha(end-i+1) = 0;
    end

    dp = zeros(size(p));
    for i = 1:size(p,1)
        % output_dp is total rate of output and input_dp is total rate of
        % input:
        output_dp = p(i) * beta(i);
        input_dp  = alpha(i);
        for j=max(1,i-(2*(fix(nuc_len/2)))):min(size(p,1),i+(2*(fix(nuc_len/2))))
            input_dp = input_dp * (1-p(j));
        end
        dp(i) = input_dp - output_dp;
    end

end

