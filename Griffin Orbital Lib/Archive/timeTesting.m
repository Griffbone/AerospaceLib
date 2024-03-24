% Frame bias matrix (from ICRS to mean equator/equinox J2000) Kaplan 2005,
% eq 3.3
da0 = -14.6/206264806.247;
xi0 = -16.6170/206264806.247;
eta0 = -6.8192/206264806.247;

B = [1 da0 -xi0; -da0, 1, -eta0; xi0, eta0, 1];

% Precession matrix
T = 0.04;
eps0 = 84381.406;

psi = ((((-0.0000000951*T + ...
	0.000132851)*T + ...
	-0.00114045)*T + ...
	-1.0790069)*T + ...
	5038.481507)*T;

w = ((((0.0000003337*T + ...
	-0.000000467)*T + ...
	-0.00772503)*T + ...
	0.0512623)*T + ...
	-0.025754)*T + eps0;

chi = ((((-0.0000000560*T + ...
	0.000170663)*T + ...
	-0.00121197)*T + ...
	-2.3814292)*T + ...
	10.556403)*T;

psi = psi*pi/(180*3600);
w = w*pi/(180*3600);
chi = chi*pi/(180*3600);

P = R3(chi)*R1(-w)*R1(eps0*pi/(180*3600))