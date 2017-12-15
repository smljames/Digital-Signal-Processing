function [X, Y] = FFTDITnbr(X, Y, N, M, KIND)
    if KIND == -1
        Y = Y*-1;
    end
    
    for m = 1:M
        L = 2^m; % L-DFT
        L2 = L/2; 
        for k = 1:L2 % k-th butterfly in L-DFT 
            Wk.r = cos(2*pi*(k-1)/L); % real part of twiddle factor in L-DFT
            Wk.i = -sin(2*pi*(k-1)/L); % imaginay part of twiddle factor in L-DFT
            
            % --- butterfly computation ---
            for i = k:L:N % top index of butterfly
                j = i + L2; % botton index of butterfly
                
                % --- multiply with twiddle factor ---
                tmp1 = X(j)*Wk.r - Y(j)*Wk.i; % real part
                tmp2 = X(j)*Wk.i + Y(j)*Wk.r; % imaginay part
                X_tmp = X(i); % previous X
                Y_tmp = Y(i); % previous Y
                X(i) = X_tmp + tmp1;
                Y(i) = Y_tmp + tmp2;
                X(j) = X_tmp - tmp1;
                Y(j) = Y_tmp - tmp2;
            end
        end
    end
    if KIND == -1
        X = X/N;
        Y = (Y*-1)/N;
    end
end