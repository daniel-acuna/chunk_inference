function chunks = create_chunks_nospace(varargin)
% Create chunks without space between them

p = inputParser;
addParamValue(p, 'n_seqlen', 10);
parse(p, varargin{:});
n_seqlen = p.Results.n_seqlen;
n_chunks = 2^n_seqlen;
chunks = zeros(n_chunks - 1, n_seqlen);
for i = 1:(n_chunks - 1)
    chunks(i, :) = bitand(i-1, 2.^(n_seqlen-1:-1:0)) > 0;
end

% Transform to numbers
chunks = cumsum(chunks, 2);

% Remove chunks that doesn't start with a number
chunks = chunks(chunks(:, 1) == 1, :);