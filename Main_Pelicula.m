% Tecnicas de Caracetrizacion, metodo de swanepoel aplicado a la
% determinacion de constantes ópticas de peliculas. 
% Juan Pablo Cruz C., Diego A. Heredia, Juan J. Medina

clear 
clc 
close all 
% se leen las funciones y se guardan en un objeto
Film = PeliculasDelgadas;

% se leen los datos de transmitancia en fuención de la longitud de onda 

%ruta donde se encuentran los datos
path = "C:\Users\Juan Pablo Nicolás\Documents\FISICA\FISICA CARRERA\10 DECIMO SEMESTRE\Tecnicas de caracterizacion\Modulo 2\Proyecto_Propiedades_Opticas\Datos.xlsx";
Data = readmatrix(path);

x=Data(1:end,1); % longitud de onda
y=Data(1:end,2); % transmitancia pelicula
ys=Data(1:end,3); % transminancia sustrato


[pksM,locsM] = findpeaks(y); % metodo para hallar los picos máximos
[pksm,locsm] = findpeaks(-y); % metodo para hallar los minimos
locsM(1)=locsM(1)-3; % correcciones a la envolvente
locsM(2)=locsM(2)-1;
locsm(1)=locsm(1)+3;


londc = [x(locsM);x(locsm)]; %longitud de onda con valores críticos

xnM = [x(locsM) ; 515.269];
ynM = [y(locsM) ; 0.117281];
TM = interp1(xnM,ynM,londc,'pchip','extrap'); % interpola los valores maximos
xnm = [x(locsm) ; 528];
ynm = [y(locsm) ; 0.172804];
Tm = interp1(xnm,ynm,londc,'pchip','extrap'); % interpola los valores minimos
s = Film.n_subst([ys(locsM);ys(locsm)]);% se calcula el indice de refracción del sustrato


Tabla0 = sortrows(table(londc,s,TM,Tm),1); % se ordenan los valores de menor a mayor longitud de onda

londc = Tabla0.londc;
TM = Tabla0.TM;
Tm = Tabla0.Tm;
s = Tabla0.s;

%plot(x,y,'*',londc,TM,'-',londc,Tm,'-');
 
n = Film.n_baja_media(s,TM,Tm); % se calcula el indice de refracción de la película

Tabla1 = sortrows(table(londc,Tm,TM,n),1);

daux = Film.espesor(londc,n); % primer calculo de los espesores
m = Film.ordenm(1237,n,londc);
df = Film.mespesor(m,londc,n);
d = Film.average(df(3:end));

%% hallando el valor de x

x_TM = Film.xTM(n,s,TM); % metodo 1

x_Tm = Film.xTm(n,s,Tm); % metodo 2

x_Ti = Film.xTi(n,s,Tm,TM); % metodo 3


alpha = Film.alpha_x(x_Ti,d);

table(londc,x_TM,x_Tm,x_Ti)

E = Film.w_to_energy(londc);

Film.tauc(E,alpha)
Film.tauc_ind(E,alpha)
k = Film.extincion(londc,alpha);

table(londc,alpha,k);

[Tteo,w] =Film.T_Teo(londc,n,s,k,d,x_Ti);
plot(x,y,'.',w,Tteo,'-',LineWidth=1.2)
xlabel('Longitud de onda, \lambda, [nm]');
ylabel('Transmitancia')
legend('Datos Experimentales','Modelo Teórico','Location','southeast')
ax = gca;
exportgraphics(ax,'modTeoExp.jpg','Resolution',300)

%{

hold on
plot(londc,Tteo)
plot(x,y)
hold off

figure 
plot(E,alpha,'*')
xlabel('Longitud de onda, \lambda, [nm]');
ylabel('\alpha, [1/nm]');
title('Coeficiente de Absorción')



figure 
plot(londc(2:end),n(2:end),'ro')
xlabel('Longitud de onda, \lambda, [nm]');
ylabel('n');
title('Índice de refracción')

figure 
plot(londc,k,'go')
xlabel('Longitud de onda, \lambda, [nm]');
ylabel('k');
title('Coeficiente de extinción')
ax = gca;
exportgraphics(ax,'extincionk.jpg','Resolution',300)

figure
inw = 1./londc(2:end).^2;
p = polyfit(inw,n(2:end),2);
y2 = polyval(p,inw);
plot(londc,n,'ro',londc(2:end),y2,'b-')
s=sprintf('n(x) = %.1f 1/x^4 + %.1f 1/x^2 + %.1f',p(1),p(2),p(3));
text(660,3,s,Color='b')
xlabel('Longitud de onda, \lambda, [nm]');
ylabel('n');
title('Índice de refracción')
legend({'n','Ajuste'})
ax = gca;
exportgraphics(ax,'indexn.jpg','Resolution',300)

%plot(londc,k,'*')
%plot(londc,n,'*')
%}
