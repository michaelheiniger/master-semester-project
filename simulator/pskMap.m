function [ C ] = pskMap( M )
% PSKMAP Creates constellation for Phase Shift Keying modulation
% C = MY_PSKMAP(M) outputs a 1×M vector with the complex symbols
% of the PSK constellation of alphabet size M, where M is an integer power of 2.

% Check if M is a power of 2 larger than 1
if (not(isNumberInteger(log2(M))) || M <= 1)
    error('M must be a power of 2 greater than 1');
end

C = exp(1i*2*pi*(0:M-1)/M);

end

