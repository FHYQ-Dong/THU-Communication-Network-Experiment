clc; clear; close all;

rng(093093);

x = zeros(1, 35);
x_interleave = interleaver(5, 7, x);
y = burst_error(x_interleave, 5);
y_deinterleave = deinterleaver(5, 7, y);

error_interleave = x_interleave ~= y;
error_deinterleave = x ~= y_deinterleave;

figure;
subplot(2,1,1);
stem(error_interleave, 'filled');
title('交织后突发错误位置');
xlabel('索引');
ylabel('错误标志');
grid on;
subplot(2,1,2);
stem(error_deinterleave, 'filled');
title('解交织后错误位置');
xlabel('索引');
ylabel('错误标志');
grid on;
