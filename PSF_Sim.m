%% PSF Simulation
% Movassagh 2024/Tuning of Freq Loop
% ------------------------------------------------------------------------
clc;
clear;
close all;
warning off;
%% ------------------------------------------------------------------------
% Read data from csv file (Train & Test)
Turbine = "Fars";
if Turbine == "Genaveh"
    data = readtable("data_genaveh/data-2024-06-11--2024-06-12.csv");
    data = data(double(strcmp(data.MBY10DU050XT03,"True"))==1,:);

elseif Turbine == "Fars"
    unit0 = "21";
    unit = "x"+ unit0;
    load data2023_02.mat
    data = data(1:10000,:);
    data = data(data.x21MBY10DU050_XT03==1,:);
end
% ----------------------------------------
% Change the header of data based on KKS and description
if Turbine == "Genaveh"
    KKS_list = [
        "MBY10DU060XT05";
        "MBY10DU050XT03";
        "MBY10DU060XT04";
        "MBY10DU060XT03";
        "MBY10DU060XT07";
        "MBY10CE901ZQ01";
        "MBY10CS901ZQ01";
        "MBY10DU050XQ01";
        "MBA11DG002XQ02";
        "MBY10FT010XQ01";
        "MBP20EU006ZV01";
        "MBP20EU007ZV01";
        "MBN20EU006ZV01";
        "MBN20EU007ZV02"];
    VariableNames = [
        "HLGIE";
        "PSFE";
        "PELRIE";
        "NTRIE";
        "ATKRIE";
        "PEL";
        "NT";
        "PSF";
        "HVL";
        "ATK";
        "NG_pre";
        "NG_diff";
        "oil_pre";
        "oil_diff"] ;

elseif Turbine == "Fars"
    KKS_list = [
        "MBY10CS901_ZQ01";
        "MBY10DS010_XQ02";
        "MBP20EU006_ZV01";
        "MBP20EU007_ZV01";
        "MBN20EU006_ZV01";
        "MBN20EU007_ZV02";
        "MBN53AA151_XQ01";
        "MBY10DU060_XQ04";
        "MBA11DG002_XQ02";
        "MBA22CT102B_XQ01";
        "MBA22CT103B_XQ01";
        "MBA22CT108B_XQ01";
        "MBA22CT104B_XQ01";
        "MBA22CT106B_XQ01";
        "MBA22CT107B_XQ01";
        "MBY10CE901_ZQ01";
        "MBY10DE010_XQ01";
        "MBY10DE010_XQ02";
        "MBY10DT010_XG01";
        "MBY10DT010_XQ07";
        "MBY10DU050_XT03";
        "MBY10DU050_XQ01";
        "MBY10DU060_XQ01";
        "MBY10DU060_XQ02";
        "MBP13AA151_XQ01";
        "MBY10DU060_XQ03";
        "MBY10DU060_XT03";
        "MBY10DU060_XT04";
        "MBY10DU060_XT05";
        "MBY10DU060_XT07";
        "MBY10FT010_XQ01";
        "MBL10CT101_XQ01";
        "MBY10FT020_XQ01";
        "MBY10FT020_XQ02";
        "MBY10FT020_XQ03";
        "MBY10FT020_XQ04";
        "MBY10FT020_XQ05";
        "MBY10FT020_XQ06";
        "MBA12CT901_ZQ01";
        "MBA11CT901_ZQ01";
        "SCADA01DE002_YB04";
        "MBY10DT010_XT01";
        "MBY10DU060_XT08";
        "MBY10DU060_XT09";
        "MBY10DT040_XQ01";
        "MBL10CM101_XQ01";
        % "MBM00EU002_XQ01";
        % "MBM00EU002_XQ02";
        % "MBM11CP103_XQ01";
        % "MBM21CP103_XQ01";
        % "MBX03CP101_XQ01";
        "MBA11CP101_XQ01";
        "MBA12CP101_XQ01";
        "MBA12CP102_XQ01";
        % "MBP15AA151_XQ01";
        % "MBN10FQ101_XQ01";
        % "MBN11CT101A_XQ01";
        % "EKG60EU101_XQ05";
        "EKG60EU101_XQ21";
        "EKG60EU101_XQ11"];

    VariableNames = [
        "NT";
        "NSV";
        "NG_pre";
        "NG_diff";
        "oil_pre";
        "oil_diff";
        "HOD";
        "HSOD";
        "HVL";
        "TET1";
        "TET2";
        "TET3";
        "TET4";
        "TET5";
        "TET6";
        "PEL";
        "PS";
        "PSV";
        "IGVplus";
        "GLT";
        "PSFE";
        "PSF";
        "YMIN";
        "JGES";
        "HG";
        "HSG";
        "NTRIE";
        "PELRIE";
        "HLGIE";
        "ATKRIE";
        "ATK";
        "TEMP_amb";
        "T_DIFF1";
        "T_DIFF2";
        "T_DIFF3";
        "T_DIFF4";
        "T_DIFF5";
        "T_DIFF6";
        "CTD";
        "CTIF";
        "SCADA_START_UP";
        "T_BLOAD";
        "VPGRIE";
        "VPVRIE";
        "YLTR";
        "HUM_amb";
        % "PRES_R1";
        % "PRES_R2";
        % "DP_CC_L";
        % "DP_CC_R"
        % 'P_HYD_FEED';
        'P_UST_COMP';
        'CPD1';
        'CPD2';
        % 'PG_VP';
        % 'FO_FEED_FLOW';
        % 'FO_temp';
        % 'FG_flow_m2';
        'FG_Temp';
        'FG_Press'] ;
    KKS_list = unit + KKS_list;
end


for i =1:length(KKS_list)
    KKS = KKS_list(i);
    KKS = KKS{1,1};
    tag = VariableNames(i);
    tag = tag{1,1};
    tag_test = tag + "_test";
    data.Properties.VariableNames = ...
        regexprep(data.Properties.VariableNames,...
        KKS, tag);
end

%% ------------------------------------------------------------------------
% First Step (call data based on variable and preprocess)

for i = 1:length(KKS_list)
    tag = VariableNames(i);
    tag = tag{1,1};
    assignin('base',tag,  fillmissing(table2array(data(:,tag)),'previous') )
    assignin('base',tag,  fillmissing(table2array(data(:,tag)),'next') )
end
if Turbine == "Genaveh"
    HLGIE = double(strcmp(HLGIE,'True'));
    PSFE = double(strcmp(PSFE,'True'));
    PELRIE = double(strcmp(PELRIE,'True'));
    NTRIE = double(strcmp(NTRIE,'True'));
    ATKRIE = double(strcmp(ATKRIE,'True'));
    NG_pre = double(strcmp(NG_pre,'True'));
    NG_diff = double(strcmp(NG_diff,'True'));
    oil_pre = double(strcmp(oil_pre,'True'));
    oil_diff = double(strcmp(oil_diff,'True'));
end
FREQUENCY = PSFE & PELRIE;
NT = NT/60;
OIL = double(oil_pre | oil_diff);
GAS = double(NG_diff | NG_pre);
t_final = length(PEL);

GLT = zeros(t_final,1);
GLT(GAS==1) = 531;
GLT(OIL==1) = 528;

WVL = GLT - 18;
VPVRIE = zeros(t_final,1);
LTRE = zeros(t_final,1);
LTRE(HLGIE~=1) = 1;

%% -------------------------------------------------------------------------
% Constants
if Turbine == "Genaveh"
    Nominal_load = 162;
    Nominal_speed = 50;
    STATNR = 0.04;
    TOTPF = 0.05;
    STAPLL=0.03;
    STAPLU = 0.08;
    STAPOM = 0.04;
    UVL25 = 0.06;
    UVL25H = 0.04;
    UVL26 = 0.02;
    UVL26H = 0.04;
    VLGMIN = 0;
    GNF2 = 38.6;
    GNF3 = 47;
    GWLF = 49.75;
    GWHF = 50.25;
    STATLF = 0.05;
    STATHF = 0.05;
    TA=1;
    KDNPF0 = 6;
    CALDN = 0.6;
    KDN = 5;
elseif Turbine == "Fars"
    Nominal_load = 162;
    Nominal_speed = 50;
    STATNR = 0.04;
    TOTPF = 0.05;
    STAPLL=0.03;
    STAPLU = 0.08;
    STAPOM = 0.04;
    UVL25 = 0.06;
    UVL25H = 0.04;
    UVL26 = 0.02;
    UVL26H = 0.04;
    VLGMIN = 0;
    GNF2 = 38.6;
    GNF3 = 47;
    GWLF = 49.75;
    GWHF = 50.25;
    STATLF = 0.04;
    STATHF = 0.04;
    TA=1;
    KDNPF0 = 6;
    CALDN = 0.6;
    KDN = 5;
end


%% ------------------------------------------------------------------------
% virtual NSV
% Notice that, this code can not make NSV when NTRIE is active. it is sth between 51.5_47.5
virtual_NSV = nan(1, t_final);
virtual_NSV(HLGIE == 1 & NT <= (2850/60)) = 50.05;
virtual_NSV(PEL >= 0.1*Nominal_load) = Nominal_speed;
% after 60s that PEL>10*Nominal PEL NSV is equal to 50
index_true_edge = find_index_true_edge(PEL >= 0.1*Nominal_load);
index_true_edge_plus60s = index_true_edge + 60;
for i = 1:length(index_true_edge)
    virtual_NSV(index_true_edge(i):index_true_edge_plus60s(i)) = nan;
end

% other signals
% VLOGNS in IGV+ is equal to 1.03, other mode is equal to 2.82, in runup ...
VLOGNS = ones(1, t_final) * 2.82;
% VLOGNS(IGV_plus == 1) = 1.03;

% VLOB
first_input = zeros(1, t_final);
second_input = zeros(1, t_final);
VLOB = zeros(1, t_final);
for i = 2:t_final
    [~, ~, first_input(i)] = LVM(NT(i), GNF3, 0, 0.001, first_input(i-1));
    [second_input(i), ~, ~] = LVM(NT(i), GNF2, 0, 0.001, second_input(i-1));
    VLOB(i) = first_input(i) && second_input(i);
end
LTRNF = ~(LTRE & ~VLOB);

% it is equal one for BUEWIE=1 or SSNOX=1
HVLABG = zeros(1, t_final);
% it is equal one when burner changeover happens or .....
STOP = zeros(1, t_final);
BLPSF = LTRNF | VPVRIE | HVLABG | STOP;

PFGLB = PELRIE | NTRIE;
LB = ones(1, t_final);
LB(NTRIE==1) = 0;


SI = zeros(1, t_final);
SI1 = zeros(1, t_final);


FGPSF = GSG0F(ATK, WVL-10, WVL-20);
% it is not equal
VLOGN5 = VLOGNS;

% Initialize variables
PS60 = virtual_NSV - NT;
PS50 = dead_line_block(PS60, TOTPF);

switch_output = ones(t_final, 1) * -1;
switch_output(PSFE == 1 & (PELRIE == 1 | ATKRIE == 1)) = 1;
PSF_OFF = zeros(t_final, 1);
PSF_OFF_Integrator_Lower = zeros(t_final, 1);

for i = 2:t_final
    [PSF_OFF(i), PSF_OFF_Integrator_Lower(i), ~] = integrator(switch_output(i), 1, 0, 30, PSF_OFF(i-1),TA);
end

PSFS = PSF_OFF .* PS50;
PS120 = PSFS / Nominal_speed;
STATPF = STAPOM;
PS130 = PS120 / STATPF;
[PS140, ~, ~] = limiter(PS130, 32/Nominal_load, -32/Nominal_load);

PS160 = GWLF - NT;
PS190 = GWHF - NT;
PS170 = zeros(t_final, 1);
PS170(PS160 > 0) = PS160(PS160 > 0);
PS200 = zeros(t_final, 1);
PS200(PS190 < 0) = PS160(PS190 < 0);
PS175 = PS170 / Nominal_speed;
PS205 = PS200 / Nominal_speed;
PS180 = PS175 / STATLF;
PS210 = PS205 / STATHF;

PS265 = max([PS140, PS180, zeros(t_final, 1)], [], 2);
PS275 = min([PS140, PS210, zeros(t_final, 1)], [], 2);
PS281 = PS275 + PS265 + 1;
PS204 = ones(t_final, 1);

C2 = ones(t_final, 1);
C3 = LTRE & PFGLB;

CF = FGPSF & ~BLPSF;
C4 = SI & LB & ~FGPSF;
SNT = SI1 | ~LB;

upper_limit = ones(t_final, 1);
lower_limit = ones(t_final, 1);

for i = 2:t_final
    if C4(i) == 0
        set_value = PS281(i);
    else
        set_value = 1;
    end

    if C2(i) == 0
        upper_limit(i) = max(1, PS204(i-1));
        lower_limit(i) = min(1, PS204(i-1));
    else
        upper_limit(i) = 3400000000000;
        lower_limit(i) = -3400000000000;
    end

    ramp_down_time = 1/0.003;
    if C3(i) == 0
        ramp_up_time = 1/0.003;
    else
        ramp_up_time = 1/0.0004;
    end

    PS204(i) = RGE(PS281(i), set_value, ramp_up_time, ramp_down_time, false, false, CF(i), SNT(i), PS204(i-1),TA);
end

Switch_val2 = (PS204 - 1) * Nominal_load;
Switch_val1 = (virtual_NSV - NT) * Nominal_load / (Nominal_speed * STATNR);
PSF_V = zeros(t_final, 1);
PSF_V(LB == 0) = Switch_val1(LB == 0);
PSF_V(LB == 1) = Switch_val2(LB == 1);

% DN
SI40 = zeros(t_final, 1);
for i = 2:t_final
    [SI40(i), ~, ~] = integrator(-SI40(i-1), 3.4*10^38, 3.4*10^-38, 5, 0,TA);
end
PSFD = SI40 + PS140 * KDNPF0 * (STATPF/KDN);
PS230 = min(PS205,PSFD);
PS220 = max (PS175 , PSFD);
PS240 = min(PS230,0) + max(PS220,0);
PS241 = PS240 + 1;

upper_limit = ones(t_final, 1);
lower_limit = ones(t_final, 1);
PS113 = ones(t_final, 1);
DNRD = 0.003 * CALDN;
DNRU = DNRD;
for i = 2:t_final
    if C4(i) == 0
        set_value = PS241(i);
    else
        set_value = 1;
    end

    if C2(i) == 0
        upper_limit(i) = max(1, PS113(i-1));
        lower_limit(i) = min(1, PS113(i-1));
    else
        upper_limit(i) = 3400000000000;
        lower_limit(i) = -3400000000000;
    end

    ramp_down_time = KDN/DNRD;
    ramp_up_time = KDN/DNRU;


    PS113(i) = RGE(PS281(i), set_value, ramp_up_time, ramp_down_time, false, false, CF(i), SNT(i), PS113(i-1),TA);
end
DN = PS113 - 1;


%% ------------------------------------------------------------------------
% Plotting
if Turbine =="Genaveh"
    figure()
    plot(PSF_V);
    hold on;
    plot(PSF);
    hold off;
    legend('PSF-V', 'PSF');
    title('Comparison of PSF-V, PSF');
    xlabel('Time');
    ylabel('Value');
elseif Turbine =="Fars"
    figure()
    plot(PSF_V);
    title('PSF-V');
    xlabel('Time');
    ylabel('Value');
end

figure()
plot(DN);
title('DN-V');
xlabel('Time');
ylabel('Value');

%% ------------------------------------------------------------------------
% Create output table
output_table = table();
output_table.PSF = PSF;
output_table.PSF_virtual = PSF_V;
output_table.PSFV_PS281 = (PS281 - 1) * Nominal_load;
output_table.Error = PSF - ((PS281 - 1) * Nominal_load);

%% ------------------------------------------------------------------------
% Write to CSV
writetable(output_table, 'output.csv');

%% ------------------------------------------------------------------------
% function

function start_indices = find_index_true_edge(boolean_series)
% Finds the starting indices of all contiguous "True" segments in a boolean series.
augmented_series = [1; double(boolean_series)]; % Convert boolean to double and prepend 1
start_indices = find(diff(augmented_series) == 1) - 1; % Find indices and adjust for MATLAB indexing
end
function output = dead_line_block(input, Threshold)
output = zeros(size(input));
output(input <= -Threshold) = input(input <= -Threshold) + Threshold;
output(input >= Threshold) = input(input >= Threshold) - Threshold;
end
function [output, QL, QU] = integrator(input, upper_lim, lower_lim, integral_time, last_output,TA)
output = last_output + (TA * input / integral_time);
if output >= upper_lim
    QU = 1;
    QL = 0;
    output = upper_lim;
elseif output <= lower_lim
    QU = 0;
    QL = 1;
    output = lower_lim;
else
    QU = 0;
    QL = 0;
end
end
function [output, QU, QL] = limiter(input, upper_lim, lower_lim)
output = input;
t_final = length(input);
QU = zeros(t_final, 1);
QL = zeros(t_final, 1);
output(input > upper_lim) = upper_lim;
output(input < lower_lim) = lower_lim;
QU(input > upper_lim) = 1;
QL(input < lower_lim) = 1;
end
function output = RGE(input, set_value, ramp_up_time, ramp_down_time, CU, CD, CF, s, last_value,TA)
if s
    output = set_value;
elseif ~CF
    output = last_value;
elseif ~s && CF && ~CU && ~CD
    if input > last_value && last_value >= 0
        output = last_value + (TA / ramp_up_time);
    elseif input < last_value && last_value <= 0
        output = last_value - (TA / ramp_up_time);
    elseif input < last_value && last_value >= 0
        output = last_value - (TA / ramp_down_time);
    elseif input > last_value && last_value <= 0
        output = last_value + (TA / ramp_down_time);
    else
        output = last_value;
    end
else
    output = last_value;
end
end
function [QU, QE, QL] = Comparitor(x1, x2)
% Ensure x1 and x2 are column vectors
x1 = x1(:);
x2 = x2(:);

% Get the length of the input vectors
t_final = length(x1);

% Initialize logical arrays
QU = false(t_final, 1);
QE = false(t_final, 1);
QL = false(t_final, 1);

% Perform comparisons
QU(x1 > x2) = true;
QE(x1 == x2) = true;
QL(x1 < x2) = true;
end
function Q = RSS(s, r, last_value)
% RS flip flop, s dominant
if ~s && ~r
    Q = last_value;
elseif ~s && r
    Q = false;
elseif s && ~r
    Q = true;
elseif s && r
    Q = true;
end
end
function QU = GSG0F(input, upper_lim, lower_lim)
% monitor the Input Value X with a Hysteresis. In a redundant
% configuration it is possible to synchronize the output signals.

% Assuming Comparitor and RSS functions are defined elsewhere

[GSG10, ~, ~] = Comparitor(input, upper_lim);
[~, ~, GSG20] = Comparitor(input, lower_lim);

% Assuming t_final is defined elsewhere in the code
GSG30 = false(size(input));

for i = 2:length(input)
    GSG30(i) = RSS(GSG10(i), GSG20(i), GSG30(i-1));
end

QU = GSG30;
end
function [QU, QM, QL] = LVM(input, average_value, interval_limit, hysteresis, last_input)
% GSG0F(input:pd.Series,upper_lim:pd.Series,lower_lim:pd.Series)

if input < average_value + interval_limit && input > average_value + interval_limit - hysteresis
    if last_input < input
        QU = false;
        QM = true;
        QL = false;
    else
        QU = true;
        QM = false;
        QL = false;
    end
elseif input < average_value - interval_limit && input > average_value - interval_limit + hysteresis
    if last_input < input
        QU = false;
        QM = false;
        QL = true;
    else
        QU = false;
        QM = true;
        QL = false;
    end
else
    QU = false;
    QM = false;
    QL = false;
end
end


