function S = sketching_mat(s, n, kind)
%SKETCHING_MAT Construct the sketch used by one solver call.

switch kind
    case 'sparse sign'
        S = sparse_sign(s, n);
    case 'exact'
        S = speye(n);
end

end
