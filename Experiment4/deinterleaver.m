function y_deinterleave = deinterleaver(row, column, y)
    y = reshape(y, row, column);
    y = y.';
    y_deinterleave = y(:).';
end
