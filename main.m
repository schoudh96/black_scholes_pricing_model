%This is the main script
%We are looking at the simulated option prices and corresponding premium of
%Asian Paints Stock Option
%Following is the historical option price data ref: 
%https://fnoanalysis.com/oi/option_chain_hist.php?symbol=ASIANPAINT&cmb_cnd_symbol=ASIANPAINT&CMB_EXPIRY_DT=2020-06-25&CMB_CND_DT=2020-05-20
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%We are simulating call option price for different strike prices and
%comparing it to the premium or running charge of buying the option
S0 = 1500; %stock price on 20th May
S_maturity = 1688.85; %stock price on 25th June
strike_prices = [1500;1520;1540;1560;1580;1600;1620;1640;1700;1760;1800];
payoff = S_maturity - strike_prices;
call_calc = zeros(11,1);
call_exact = zeros(11,1);
abs_error = zeros(11,1);
premium = [74.85;65.00;60.00;51.40;57.20;36.65;35;37.40;19.80;15.00;11.10];
for i = 1:11
    call_calc(i) = Black_Scholes_Script(S0,strike_prices(i));
    call_exact(i) = exact_call_price(S0,strike_prices(i));
    abs_error(i) = abs(call_exact(i) - call_calc(i));
end
analysis_table = table(strike_prices,premium, call_exact, call_calc, abs_error,payoff);
disp(analysis_table);
analysis_table.Variables = round(analysis_table.Variables,5);
writetable(analysis_table,'analysis_table.csv','Delimiter',',');