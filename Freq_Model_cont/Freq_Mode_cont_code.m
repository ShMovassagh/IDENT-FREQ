%% load Data

load 'D:\khani\Documents\3-Digitalization\Gas Turbine\ident\laptop\GT21\Speed_load_Controller\GT21_2023_F_PSF_V1'
n=[535000  617000;    1900000 1975000;    646000  773500;    1179000 1420000;    197400  309200;    606000  652000;    1381900 1412500;    184000  207000;    520000  546400;    928000  1200000; 4000     210000;    1480000 1540000;
    ];
for i=1:12
    n11=n(i,1);
    n12=n(i,2);
    dd{i}=data{i}(n11:n12,:);
    dd{i}.NSTERN=nstern_f(dd{i}.CTIF_av_in,dd{i}.P_UST_COMP,dd{i}.NT,3000);
    du{i}=[dd{i}.YMIN  dd{i}.CPD1-1.01325 dd{i}.CTD_av_out dd{i}.CTD_av_out.^2 dd{i}.NSTERN];
    dy{i}=dd{i}.PEL2;
end

%Preparing state space representation of identified model
load 'D:\khani\Documents\3-Digitalization\Gas Turbine\ident\laptop\GT21\Models_Final\U27_month11&5'
Model=sysP1D;
m=tf(Model);
m2=ss(m);
X(:,1) = zeros(1,size((m2.A),2));

%% Parameters
XDL=-1;
XDU=2;


PN=162;
NSV=3000/60;
NN=3000/60;
Droop=0.05;
KDN=5;

% Preparing inputs & outputs (month #j)
j = 11;
U_r = du{j}(:,2:5)';
YMIN = du{j}(:,1);
PEL_act = dy{j};
PSWN = dd{j}.PSV/PN;
PELN = dd{j}.PEL/PN;
PSFN = dd{j}.PSF_V/PN;
DN = dd{j}.DN_V;

n=size(PSWN,1);

% Controller Parameter
kp = 0.1;
TN = 1;
alfa = 0.5191;


YI(1)=0;
%% Control Loop
for i = 1:size(PSWN,1)
    e(i) = PSWN(i)-PELN(i);

    %Saturation #1
    if e(i) > XDU
        e(i) = XDU;
    elseif e(i) < XDL
        e(i) = XDL;
    end

    XD(i) = e(i)+PSFN(i);
    Yp(i) = kp*XD(i);

    %Integral
    YI(i+1) = YI(i)+(1/TN)*Yp(i);

    FF(i) = KDN * DN(i) + alfa * PSWN(i);

    % Saturation #2 Upper limit
    YU(i) = min(1,(8*XD(i)+YMIN(i))) - FF(i) -Yp(i);
    % Saturation #2 Lower limit
    YL(i) = -1 - Yp(i);

    %Saturation #2
    if YI(i+1)> YU(i)
        YI(i+1) = YU(i);
    elseif YI(i)< YL(i)
        YI(i+1) = YL(i);
    end

    Y(i)=YI(i+1)+Yp(i);

    YMIN(i) = 100* (max(Y(i)+Yp(i),-1)+ FF(i));

    % Model
    U(:,i)=[YMIN(i);U_r(:,i)];
    X(:,i+1)=(m2.A+eye(size(m2.A)))*X(:,i)+m2.B*U(:,i);
    PEL(i)=m2.C*X(:,i+1);%+m4.D*U(:,i);

    PELN(i) = PEL(i)/PN;
end

figure();plot(PEL);hold on;plot(PEL_act);


%% check Model
% for i=1:n
%     U(:,i)=du{j}(i,:)';
%     X(:,i+1)=(m2.A+eye(size(m2.A)))*X(:,i)+m2.B*U(:,i);
%     PEL(i)=m2.C*X(:,i+1);%+m4.D*U(:,i);
% end
% figure();plot(PEL);hold on;plot(dy{j})

%%