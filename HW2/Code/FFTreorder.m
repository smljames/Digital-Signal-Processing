function x = FFTreorder(x)
    N = length(x);
    M = log2(N);
    if M == 1
        % return x
    else
    x_1 = x(1:2:end);
    x_2 = x(2:2:end);
    x_1 = FFTreorder(x_1);
    x_2 = FFTreorder(x_2);
    x = [x_1 x_2];
    end
end