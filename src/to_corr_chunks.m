function cor_chunks = to_corr_chunks(chunks)
% Transform chunks into correlation chunks
cor_chunks = zeros(size(chunks, 2), size(chunks, 2), size(chunks, 1));
for z = 1:size(chunks, 1)
    cur_d = 0;
    for j = 1:size(chunks, 2)
        if chunks(z, j) == 0
            cur_d = 0;
            cor_chunks(j, j, z) = 1;
        elseif (j>1 && chunks(z, j) ~= chunks(z, j-1))
            cur_d = 1;
            cor_chunks(j, j, z) = 1;
        else
            cor_chunks(max(j-cur_d, 1):j, ...
                max(j-cur_d, 1):j, z) = 1;
            cur_d = cur_d + 1;
        end
    end
end