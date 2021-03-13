close all;
clear;
clc;


tab = readtable("../radiation_patterns/phi=0.txt");
RP_single = tab{:,3};
tmp = flip(RP_single,1);
tmp = tmp(1:end-1)';

RP_db = [tmp RP_single(1:end-1)'];
%normalize RP
RP_db_max = max(RP_db);
RP_db = RP_db-RP_db_max;

%convert RP to voltage
RP_volt = 10.^(RP_db/20);

plot(RP_volt);