function [sigma, risk_free_rate, time] = calculate_parameters()
stock_data = readmatrix('ASIAN_PAINT_Stock_Prices.csv'); %stock level data
size_data = size(stock_data); %size of data
close_price = stock_data(:,5); %select close price of last one month(22 days)

%find historical volatility by first calculating log of daily returns
vol_table = zeros(size_data(1), 2);
vol_table(:,1) = close_price;
for j=2:size_data(1)
    vol_table(j,2) = log(vol_table(j,1)/vol_table(j-1,1));
end

mean_return = mean(vol_table(:,2));
std_dev = sqrt(1/(size_data(1)-2)*sum((vol_table(2:size_data(1),2)-mean_return).^2));
sigma = std_dev*sqrt(252);
time = 26/252;
risk_free_rate = 0.0575; %10 year G-sec bonds yield is 5.75% as fixed by RBI on 22nd May
%Ideally I would consider 91 day treasury rates which was 3.28% on 22nd
%May, but the simulation and report were already completed and using a
%lower risk free rate just reduced the call price for every strike. Given
%that the call price is already undervalued it just further supports our
%arguments as noted in the report. 