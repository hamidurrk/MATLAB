time = [0, 1, 2, 3, 4, 5, 6]; 
A_conc = [5, 4.093654, 3.3516, 2.744058, 2.246645, 1.839397, 1.505971];
ln_A = log(A_conc);

p = polyfit(time, ln_A, 1);
slope = p(1);
k = -p(1); 
disp(['Estimated k: ', num2str(k)]);

plot(time, ln_A, 'o-', 'LineWidth', 2, 'MarkerSize', 5); 
xlabel('Time, t (s)'); 
ylabel('ln[A] (Natural Log of Concentration)'); 
title('First-Order Reaction: ln[A] vs. Time');
grid on;

x_pos = max(time) - 3; 
y_pos = max(ln_A) - 0.1;
text(x_pos, y_pos, ['Slope = ', num2str(slope, '%.4f')], 'FontSize', 12, 'Color', 'red');
x_pos = max(time) - 3; 
y_pos = max(ln_A) - 0.2;
text(x_pos, y_pos, ['k = ', num2str(k, '%.4f')], 'FontSize', 12, 'Color', 'red')
