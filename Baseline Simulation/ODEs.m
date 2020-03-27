function [dydt] = ODEs(t,y,params)
% ODEs defines the system of ODEs describing the model

% assign parameter values
f1 = params(1);%env-IgG1 kon forward
r1 = params(2);%koff reverse
f2 = params(3);%env-IgG2
r2 = params(4);
f3 = params(5);%env-IgG3
r3 = params(6);
f4 = params(7);%env-IgG4
r4 = params(8);
f5 = params(9);%IgG1-FcR
r5 = params(10);
f6 = params(11);%IgG2-FcR
r6 = params(12);
f7 = params(13);%IgG3-FcR
r7 = params(14);
f8 = params(15);%IgG4-FcR
r8 = params(16);
g1tot    = params(17);%mM 
g2tot    = params(18);%mM 
g3tot    = params(19);%mM 
g4tot    = params(20);%mM 
etot    = params(21);%mM 
ftot    = params(22);%mM 

% assign complex concentrations 
e1 = y(1); %env+IgG1
e2 = y(2);  
e3 = y(3); 
e4 = y(4);
e11 = y(5); %env+IgG1+IgG1
e12 = y(6);
e13 = y(7); 
e14 = y(8); 
e22 = y(9); 
e23 = y(10);   
e24 = y(11);
e33 = y(12);
e34 = y(13);
e44 = y(14);
fe11 = y(15); %FcR dimer+IgG+IgG+env
fe12 = y(16);
fe13 = y(17); 
fe14 = y(18); 
fe22 = y(19); 
fe23 = y(20);
fe24 = y(21);
fe33 = y(22);
fe34 = y(23);
fe44 = y(24);

% on an off rates for each specific reaction
k1f = f1; %kon for reaction 1
k1r = r1; %koff for reaction 1
k2f = f2; %kon for reaction 2 ...
k2r = r2;
k3f = f3;
k3r = r3;
k4f = f4;
k4r = r4;
k5f = f1;
k5r = r1;
k6f = f2;
k6r = r2;
k7f = f1;
k7r = r1;
k8f = f3;
k8r = r3;
k9f = f1;
k9r = r1;
k10f = f4;
k10r = r4;
k11f = f1;
k11r = r1;
k12f = f2;
k12r = r2;
k13f = f3;
k13r = r3;
k14f = f2;
k14r = r2;
k15f = f4;
k15r = r4;
k16f = f2;
k16r = r2;
k17f = f3;
k17r = r3;
k18f = f4;
k18r = r4;
k19f = f3;
k19r = r3;
k20f = f4;
k20r = r4;
k21f = (f5+f5)/2; %kon/off for reaction 21 - 30 are the average of the two 
k21r = (r5+r5)/2; %kon/offs for the two IgGs involved in the reaction
k22f = (f5+f6)/2;
k22r = (r5+r6)/2;
k23f = (f5+f7)/2;
k23r = (r5+r7)/2;
k24f = (f5+f8)/2;
k24r = (r5+r8)/2;
k25f = (f6+f6)/2; 
k25r = (r6+r6)/2;
k26f = (f6+f7)/2;
k26r = (r6+r7)/2;
k27f = (f6+f8)/2;
k27r = (r6+r8)/2;
k28f = (f7+f7)/2;
k28r = (r7+r7)/2;
k29f = (f7+f8)/2; 
k29r = (r7+r8)/2;
k30f = (f8+f8)/2; 
k30r = (r8+r8)/2;

% Conservation equations
g1 = g1tot-(e1+2*e11+e12+e13+e14+2*fe11+fe12+fe13+fe14);
g2 = g2tot-(e2+2*e22+e12+e23+e24+2*fe22+fe12+fe23+fe24);
g3 = g3tot-(e3+2*e33+e23+e13+e34+2*fe33+fe23+fe13+fe34);
g4 = g4tot-(e4+2*e44+e24+e34+e14+2*fe44+fe24+fe34+fe14);
e = etot-(e1+e2+e3+e4+e11+e12+e13+e14+e22+e23+e24+e33+e34+e44+fe11+fe12+...
    fe13+fe14+fe22+fe23+fe24+fe33+fe34+fe44);      
f = ftot-(fe11+fe12+fe13+fe14+fe22+fe23+fe24+fe33+fe34+fe44);
                   
% Reaction rates
react1  = k1f*g1*e-k1r*e1;
react2  = k2f*g2*e-k2r*e2;
react3  = k3f*g3*e-k3r*e3;
react4  = k4f*g4*e-k4r*e4; 
react5  = k5f*g1*e1-k5r*e11;
react6  = k6f*g2*e1-k6r*e12;
react7 = k7f*g1*e2-k7r*e12;
react8 = k8f*g3*e1-k8r*e13;
react9 = k9f*g1*e3-k9r*e13;
react10 = k10f*g4*e1-k10r*e14;
react11 = k11f*g1*e4-k11r*e14;
react12 = k12f*g2*e2-k12r*e22; 
react13 = k13f*g3*e2-k13r*e23;
react14 = k14f*g2*e3-k14r*e23;
react15 = k15f*g4*e2-k15r*e24;
react16 = k16f*g2*e4-k16r*e24;
react17 = k17f*g3*e3-k17r*e33;
react18 = k18f*g4*e3-k18r*e34;
react19 = k19f*g3*e4-k19r*e34;
react20 = k20f*g4*e4-k20r*e44;
react21 = k21f*e11*f-k21r*fe11;
react22 = k22f*e12*f-k22r*fe12; 
react23 = k23f*e13*f-k23r*fe13;
react24 = k24f*e14*f-k24r*fe14;
react25 = k25f*e22*f-k25r*fe22;
react26 = k26f*e23*f-k26r*fe23;
react27 = k27f*e24*f-k27r*fe24;
react28 = k28f*e33*f-k28r*fe33;
react29 = k29f*e34*f-k29r*fe34;
react30 = k30f*e44*f-k30r*fe44;

% Differential equations;
de1 = react1-react5-react6-react8-react10;
de2 = react2-react7-react12-react13-react15;
de3 = react3-react9-react14-react17-react18;
de4 = react4-react11-react16-react19-react20;
de11 = react5-react21;
de12 = react6+react7-react22;
de13 = react8+react9-react23;
de14 = react10+react11-react24;
de22 = react12-react25;
de23 = react13+react14-react26;
de24 = react15+react16-react27;
de33 = react17-react28;
de34 = react18+react19-react29;
de44 = react20-react30;
dfe11 = react21;
dfe12 = react22;
dfe13 = react23;
dfe14 = react24;
dfe22 = react25;
dfe23 = react26;
dfe24 = react27;
dfe33 = react28;
dfe34 = react29;
dfe44 = react30;

% [uM/s] product
dydt=[de1;de2;de3;de4;de11;de12;de13;de14;de22;de23;de24;de33;de34;de44;...
    dfe11;dfe12;dfe13;dfe14;dfe22;dfe23;dfe24;dfe33;dfe34;dfe44];  
% Reassemble differential equations

% Manipulation of global variables for storing fluxes, intermediate variables
global tStep tArray varArray;
if t > tArray(tStep)            % Roughly eliminates data from rejected time steps
    tStep = tStep + 1;  
 end
tArray(tStep) = t;              
varArray(tStep,:) = [g1,g2,g3,g4,e,f];  
end