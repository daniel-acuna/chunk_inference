function res = apply(func, X, varargin)
% Apply function row-wise to X
p = inputParser;
addRequired(p, 'func', @(x)(isa(x, 'function_handle')));
addRequired(p, 'X', @isnumeric);
addOptional(p, 'dim', 2);
parse(p, func, X, varargin{:});
dim = p.Results.dim;
C = num2cell(X, dim);
res = cell2mat(cellfun(func, C, 'UniformOutput', false));