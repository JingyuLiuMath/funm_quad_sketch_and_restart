function [ell, U, T, evT] = thick_quad(U, T, evT, number_thick)
%THICK_QUAD Retain modes with the smallest absolute real parts.
% A selected real Schur block is retained in full, including conjugate pairs.

number_thick = min(number_thick, numel(evT));
[~, ind] = sort(abs(real(evT)), 'ascend');
selected = false(size(evT));
selected(ind(1:number_thick)) = true;

if isreal(T)
    j = 1;
    while j < size(T, 1)
        if T(j + 1, j) ~= 0
            selected(j:j + 1) = any(selected(j:j + 1));
            j = j + 2;
        else
            j = j + 1;
        end
    end
end

ell = sum(selected);
[U, T] = ordschur(U, T, selected);
evT = ordeig(T);

end
