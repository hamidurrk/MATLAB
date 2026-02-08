function [yearly_sales, yearly_profit, quarterly_profit] = ex3p5()
    load('business.mat', 'cost', 'price', 'sales');
    
    yearly_sales = sum(sales, 2); % 2x1 vector
    
    profit_per_quarter = (price - cost) .* sales; % 2x4 matrix
    yearly_profit = sum(profit_per_quarter, 2);   % 2x1 vector
    
    quarterly_profit = sum(profit_per_quarter, 1); % 1x4 vector
end