%% Identificare date

t = scope115(:,1);
u = scope115(:,2);
y1 = scope115(:,3);
y2 = scope115(:,4);

plot(t,[u y1]);
xlabel('Timp [s]');
ylabel('u,y1 [ V ]');
legend('u','y1');
shg;title("Datele masurate");grid;
set(gcf, 'Color', 'w');

%% Determinarea functiei de transfer de gradul 2
Mr = (y1(447)-y1(459)) / (u(443)-u(456)); 
Tr = 2*(t(459) - t(447)); 
wr = 2*pi/Tr; % frecventa de rezonanta
tita=sqrt(2-sqrt(4-4/Mr^2))/2;
wn = wr/sqrt(1-2*tita^2);
K = mean(y1)/mean(u);
H = tf(K*(wn^2),[1 2*tita*wn wn^2]);
y1c = lsim(H,u,t);

%% Spatiul starilor pentru conditi initiale nenule
A = [0 1; -wn^2 -2*tita*wn];
B = [0; K*wn^2];
C = [1 0];
D = 0;
y1css = lsim(A,B,C,D,u,t,[y1(1) 0]);
figure;
plot(t,y1,t,y1css);shg
legend('y1','y1calculat');
set(gcf, 'Color', 'w');
xlabel('Timp [s]');
ylabel('y1,y1calculat ');
title("Sistemul obtinut dupa identificare")
eMPN = norm(y1css-y1)/norm(y1css-mean(y1))
eMPN2 = norm(y1-y1css)/norm(y1-mean(y1))

%% Raspunsul in frecventa
figure;
w = logspace(2,5);
[num,den] = tfdata(H,'v');
[M,Ph] = bode(num,den,w);


% la Frecventa joasa
w1 = pi/(t(173)-t(133));
w2 = pi/(t(229)-t(201));
w3 = pi/(t(274)-t(252));
w4 = pi/(t(313)-t(294));


M1 = (y1(173)-y1(133))/(u(167)-u(128));
M2 = (y1(229)-y1(201))/(u(226)-u(201));
M3 = (y1(252)-y1(274))/(u(250)-u(271));
M4 = (y1(294)-y1(313))/(u(292)-u(310));

% La rezonanta
w5 = pi/(t(570)-t(559));
M5 = (y1(559)-y1(570))/(u(555)-u(565));

w6 = pi/(t(678)-t(669));
M6 = ((y1(669)-y1(678))/(u(664)-u(673)));


% la Frecventa Inalta
w7= pi/(t(854)-t(847));
w8 = pi/(t(929)-t(923));

M7 = (y1(847)-y1(854))/(u(842)-u(849));
M8 = (y1(923)-y1(929))/(u(918)-u(924));



omega = [w1 w2 w3 w4 wr w5 w6 w7 w8]
modul = [M1 M2 M3 M4 Mr M5 M6 M7 M8];

subplot(211);
semilogx(omega,20*log10(modul),'x');grid;
title("Caracteristica de modul");
xlabel("w(lg)");
ylabel("|M|^d^B")
hold on;set(gcf, 'Color', 'w');
plot(w,20*log10(M))

% Faza la frecventa joasa



ph1 = rad2deg((t(128)-t(133))*w1);
ph2 = rad2deg((t(226)-t(229))*w2);
ph3 = rad2deg((t(271)-t(274))*w3);
ph4 = rad2deg((t(310)-t(313))*w4);

% Faza la rezonanta
w_rez = wr;
phr = rad2deg((t(443)-t(447))*wr);
ph5 = rad2deg((t(565)-t(570))*w5);
ph6 = rad2deg((t(673)-t(678))*w6);

% Faza la frecvente inalte
ph7 = rad2deg((t(849)-t(854))*w7);
ph8 = rad2deg((t(924)-t(929))*w8);

omega_ph = [w1 w2 w3 w4 wr w5 w6 w7 w8];
faza = [ph1 ph2 ph3 ph4 phr ph5 ph6 ph7 ph8];
subplot(212);
semilogx(omega_ph,faza,'x');grid;
title("Caracteristica de faza");
xlabel("w(lg)");
ylabel("Faza in grade")
hold on;
plot(w,Ph)
set(gcf, 'Color', 'w');

%% Partea 2 - Identificarea Parametrica
%% Y1
%% ARMAX 
dt = t(2)-t(1); 
d_y = iddata(y1,u,dt);
M_y1_armax = armax(d_y,[2 2 2 1]);
num = M_y1_armax.B;
den = M_y1_armax.A;
Hz_armax = tf(num,den,dt,'Variable','z^-1');
compare(d_y,M_y1_armax)
figure;
resid(d_y,M_y1_armax,7);

%% verificare
y1c = dlsim(Hz_armax,t,u);
compare(y1c,iddata(y1))
eMPN2 = norm(-y1c+y1)/norm(y1-mean(y1))
eMPN1 = norm(y1c-y1)/norm(y1c-mean(y1))


%% Oe
dt = t(2)-t(1);
d_y = iddata(y1,u,dt);
M_y1_oe = oe(d_y,[2 2 0]);
num = M_y1_oe.B;
den = M_y1_oe.F;
Hz_oe = tf(num,den,dt,'Variable','z^-1');
compare(d_y,M_y1_oe);
figure;
resid(M_y1_oe,d_y,7)

%%
figure;
y1c = dlsim(Hz_oe,t,u);
compare(y1c,iddata(y1))





 

%% Y2
%% ARMAX
dt = t(2)-t(1);
d_y = iddata(y2,u,dt);
M_y2_armax = armax(d_y,[2 2 2 0]);
num = M_y2_armax.B;
den = M_y2_armax.A;
Hz2_armax = tf(num,den,dt,'variable','z^-1'); %% confirmare cu dlsim sau lsim
compare(d_y,M_y2_armax)
figure;
resid(d_y,M_y2_armax,'corr',7);


y2c = dlsim(Hz2_armax,t,u);
figure;
compare(y2c,iddata(y2))



%% Oe 
dt = t(2)-t(1);
d_y = iddata(y2,u,dt);
M_y2_oe = oe(d_y,[1 2 0]);
num = M_y2_oe.B;
den = M_y2_oe.F;
Hz2_oe = tf(num,den,dt,'variable','z^-1');
compare(d_y,M_y2_oe);
figure;
resid(M_y2_oe,d_y,'corr',5);

figure;
[num,den] = tfdata(Hz2_oe,'v');
y1c = dlsim(num,den,u);
compare(y1c,iddata(y2));

%%
num = [0.108 0 0];
den = [1 -1.659 0.7631];
y1c = dlsim(num,den,u);
compare(y1c,iddata(y2))




