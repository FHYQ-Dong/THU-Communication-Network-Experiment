clc; clear; close all;

rng(093);

x = 1:35;
x_interleave = interleaver(5, 7, x);
x_deinterleave = deinterleaver(5, 7, x_interleave);
figure;
subplot(3,1,1);
stem(x, 'filled');
title('原始序列');
xlabel('索引');
ylabel('数值');
grid on;
subplot(3,1,2);
stem(x_interleave, 'filled');
title('交织序列 (5x7)');
xlabel('索引');
ylabel('数值');
grid on;
subplot(3,1,3);
stem(x_deinterleave, 'filled');
title('解交织序列');
xlabel('索引');
ylabel('数值');
grid on;
