clear all; %close all; clc;

%%% This script defines the parameters that need to be defined and then
%%% runs all necessary subscripts to be able to output and plot the
%%% performance of the rocket with time. This script should make it very
%%% easy to modify input parameters to observe their effects on the output
%%% parameters so we can tweak the performance to match that which is
%%% desired.

%%% First created by Dev on 23 Feb 2018

F_range=150:100:250;
I_range=1000:100:1200;

for i=1:length(F_range)
for     j=1:length(I_range)
    
[I_total_result, t_burn_result,F_avg_result,Isp_avg, m_f_total, m_ox_total, Lp,Dia_init,Dia_end] = makeArocketLoop( I_range(j), F_range(i));

I_total_store(i,j) = I_total_result;
end
end

