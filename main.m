Ts = 0.001;

m1  = 20;
c1  = 0.1;
Ka1 = 10;
Kt1 = 1;
Kp1 = 200;
Ki1 = 0;
Kd1 = 20;

Gol1 = tf([Ka1*Kt1],[m1, c1, 0]);

COF = 60*2*pi;

[mag,phi] = bode(Gol1, COF);
phim = 60-180-phi;
phim = phim*pi/180;
alpha = (1+sin(phim))/(1-sin(phim));
T = 1/(sqrt(alpha)*COF);

Cll1num = [alpha*T, 1];
Cll1den = [T,1];

Cll1 = tf(Cll1num, Cll1den);

[mag,phi] = bode(Cll1*Gol1, COF);
K1 = 1/mag;

NLT = K1*Cll1*Gol1;

m2  = 0.2;
c2  = 0.01;
Ka2 = 10;
Kt2 = 1;
Kp2 = 0.1;
Ki2 = 0;
Kd2 = 0;

Gol2 = tf([Ka2*Kt2],[m2, c2, 0]);

COF = COF*10;

[mag,phi] = bode(Gol2, COF);
phim = 60-180-phi;
phim = phim*pi/180;
alpha = (1+sin(phim))/(1-sin(phim));
T = 1/(sqrt(alpha)*COF);

Cll2num = [alpha*T, 1];
Cll2den = [T,1];

Cll2 = tf(Cll2num, Cll2den);

[mag,phi] = bode(Cll2*Gol2, COF);
K2 = 1/mag;

NLT = K2*Cll2*Gol2;

E = 1/(1+Cll1*K1*Gol1);

Kpnm = 50;
Kdnm = 1;