function [ result ] = isNumberInteger(number)
% IsNumberInteger
%  result = isNumberInteger(number) return 1 if number is an integer and 0 else

result = (mod(number,1) == 0);

end

