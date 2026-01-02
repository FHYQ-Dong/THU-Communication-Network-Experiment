function x_interleave = interleaver(row, column, x)
    x = reshape(x, column, row);
    x = x.';
    x_interleave = x(:).';
end
