function y = burst_error(x, L)
    noise = zeros(1,size(x,2));
    error_idx = randi([1,size(x,2)-L+1]);
    noise(1,error_idx:error_idx+L-1) = (rand(1,L)<0.5);
    y = mod(x + noise, 2);
end
