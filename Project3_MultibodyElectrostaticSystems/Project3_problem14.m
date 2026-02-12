clc
clear

a = 0.1e-3;
h = 1e-3;
d = 0.5e-3;
epsilon_r = 4;
epsilon0 = 8.854e-12;
epsilon = epsilon0 * epsilon_r;
V2 = 2;
V3 = -2;

a11 = 1/(2*pi*epsilon)*log(2*h/a);
a12 = 1/(2*pi*epsilon)*log(sqrt((2*h)^2+d^2)/d);
a13 = 1/(2*pi*epsilon)*log(sqrt((2*h)^2+(2*d)^2)/(2*d));
a14 = 1/(2*pi*epsilon)*log(sqrt((2*h)^2+(3*d)^2)/(3*d));

a21 = a12;
a22= a11;
a23 = 1/(2*pi*epsilon)*log(sqrt((2*h)^2+d^2)/d);
a24 = 1/(2*pi*epsilon)*log(sqrt((2*h)^2+(2*d)^2)/(2*d));

a31 = a13;
a32 = a23;
a33=a11;
a34 = 1/(2*pi*epsilon)*log(sqrt((2*h)^2+(d)^2)/(d));

a41 = a14;
a42 = a24;
a43 = a34;
a44=a11;

a= [a11 a12 a13 a14;
    a21 a22 a23 a24;
    a31 a32 a33 a34;
    a41 a42 a43 a44];

a_q2 = a22 - a32;
a_q3 = a23 - a33;

Q_2 = V2 / (a22+a23);

V_1 = (a12-a13)*Q_2;
V_4 = (a42 - a43) * Q_2;