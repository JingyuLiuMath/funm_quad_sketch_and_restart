function rho = evalnodal(z, nodes, subdiag)

rho = ones(size(z));
for j = 1 : length(nodes)
    rho = rho * subdiag(j) ./ (z - nodes(j));
end

end
