function [X, Y] = FFTDIF(X, Y, N, M, KIND)
    if KIND == -1
        Y = Y*-1;
    end
    
    for m = M:-1:1
        L = 2^m; % L-DFT
        L2 = L/2; 
        for k = 1:L2 % k-th butterfly in L-DFT 
            Wk.r = cos(2*pi*(k-1)/L); % real part of twiddle factor in L-DFT
            Wk.i = -sin(2*pi*(k-1)/L); % imaginay part of twiddle factor in L-DFT
            
            % --- butterfly computation ---
            for i = k:L:N % top index of butterfly
                j = i + L2; % botton index of butterfly
                
                % --- multiply with twiddle factor ---
                Xi_tmp = X(i); % previous X(i)
                Yi_tmp = Y(i); % previous Y(i)
                Xj_tmp = X(j); % previous X(j)
                Yj_tmp = Y(j); % previous Y(j)
                X(i) = Xi_tmp + Xj_tmp; 
                Y(i) = Yi_tmp + Yj_tmp;
                X(j) = (Xi_tmp - Xj_tmp) * Wk.r - (Yi_tmp - Yj_tmp) * Wk.i;
                Y(j) = (Yi_tmp - Yj_tmp) * Wk.r + (Xi_tmp - Xj_tmp) * Wk.i;
            end
        end
    end
    
    if KIND == -1
        X = X/N;
        Y = (Y*-1)/N;
    end
    
    X = FFTreorder(X);
    Y = FFTreorder(Y);
end