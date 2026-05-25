function [idx, val] = middleValue(A)
%MIDDLEVALUE  Middle index and value of an array.
%
%   [IDX, VAL] = MIDDLEVALUE(A) returns the middle index IDX and the
%   corresponding value VAL from array A.
    arguments
        A {mustBeNonempty}
    end

    n = numel(A);

    if mod(n,2) == 0
        idx = n/2;
    else
        idx = (n+1)/2;
    end

    val = A(idx);
end
