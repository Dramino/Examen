%% Programa para obtener la grafica de velocidades
% FlutterVcr.m
clc
clear all
[H1i,H2i,H3i,H4i,A1i,A2i,A3i,A4i]=derAer(1,2,200);
v=linspace(1,2,200);

[archivo,direccion] = uigetfile({'*.txt'},'File Selector');

dirCom     = strcat(direccion,archivo);
archivoID  = fopen(dirCom,'r');
datos      = textscan(archivoID, '%s','Delimiter', '\n');
matriz     = datos{1};

varConst  = zeros(1,7);
revConst  = false;
revMat    = false;
numConst  = 0;
matConst  = 0;

for i = 1 : length(matriz)
  if revConst == false
    constante = strsplit(matriz{i},'=');
    if length (constante)== 2
      num                = str2double(strsplit(constante{2}));
      numConst           = numConst + 1;
      varConst(numConst) = num(2);
      if numConst == 7, revConst=true; end
    end
  elseif revMat == false && revConst==true
        mat    = strsplit(matriz{i});
        matNum = str2double(mat);
        if not(isnan(sum(matNum)))
            switch matConst
            case 0
              L        = matNum;
              matConst = matConst + 1;
            case 1
              fiZ      = matNum;
              matConst = matConst + 1;
            case 2
              fiT      = matNum;
              matConst = matConst + 1;
              revMat   = true;
            end
        end
  end
end

mz  = varConst(1);
mt  = varConst(2);
B   = varConst(3);
wz  = varConst(4);
wt  = varConst(5);
zz  = varConst(6);
zt  = varConst(7);
rho = 1.25;

%Integral de formas modales
intFiZ  = integrarS(fiZ.^2,L);
intFiT  = integrarS(fiT.^2,L);
intFiZT = integrarS(fiZ.*fiT,L);

%cambio de variable
g  = wt/wz;
Bz = rho*B^2/mz;
Bt = rho*B^4/mt;
Y  = intFiZT^2/(intFiT*intFiZ);

dominio = length(v);
GrafRe  = zeros(1,dominio);
GrafIm  = zeros(1,dominio);

%% Obtencion de raices
for i=1:dominio
    H1=H1i(i);
    H2=H2i(i);
    H3=H3i(i);
    H4=H4i(i);

    A1=A1i(i);
    A2=A2i(i);
    A3=A3i(i);
    A4=A4i(i);

  %Reales
  R0 = 1;
  R1 = 0;
  R2 = -(1+g^2+4*g*zz*zt+Bz/2*g^2*H4+Bt/2*A3);
  R3 = g*(zt*Bz*g*H1+zz*Bt*A2);
  R4 = g^2*(1+Bz/2*H4+Bt/2*A3+Bz*Bt/4*(A1*H2*Y-A2.*H1+A3*H4-A4*H3*Y));

  raices    = roots([R4 R3 R2 R1 R0]);
  r         = raices(imag(raices)==0);
  GrafRe(i) = max(r);

  %Imaginarios
  I0 = -2*(+zz*g+zt);
  I1 = -2*(1/4*(Bz*g^2*H1+Bt*A2));
  I2 = 2*((zz*(Bt/2*A3+g)+zt*g^2*(Bz/2*H4+1)));
  I3 = 2*(g^2*((Bz*Bt)/8*(H1*A3*Y-H2*A4*Y-H3*A1+H4*A2)+1/4*(Bz*H1+Bt*A2)));

  raices    = roots([I3 I2 I1 I0]);
  r         = raices(imag(raices)==0);
  GrafIm(i) = max(r);

end
%% Grafica
fig1 = figure(1);
plot    (v,GrafRe,'b-',v,GrafIm,'r-','LineWidth',1.25);
title   ('Solucion para las ecuaciones de aleteo')
xlabel  ('$\hat{V}_{cr}$','Interpreter','latex','Fontsize',11)
ylabel  ('$\hat{\omega}_r$','Interpreter','latex','Fontsize',11)
grid on
h = legend('Raices reales','Raices imaginarias');
set(h,'Interpreter','latex','Location','SouthWest')

legend('boxoff')

%% Punto de intesreccion
vcr=0;wr=0;
if GrafRe(1) - GrafIm(1) >= 0
  dif=1;
else
  dif=-1;
end
for i=1:length(GrafRe)
  if (GrafRe(i) - GrafIm(i)) * dif < 0
    [vcr,wr]=interseccion(GrafRe(i-1),GrafRe(i),GrafIm(i-1),GrafIm(i),v(i-1),v(i));
  end
end

%% Resultados
archivo= fopen ('Resultado.txt', 'w');
fprintf(archivo,'%s\r\n',datetime);
fprintf(archivo,'Velocidad critica: %6.4f\r\n',vcr);
fprintf(archivo,'Frecuencia: %6.4f\r\n\r\n',wr);
fprintf(archivo,'%6s %6s %6s\r\n','Vcr', 'Re', 'Im');
fprintf(archivo,'%6.4f %6.4f %6.4f\r\n',v, GrafRe, GrafIm);
fclose(archivo);

%% Funciones
function [x,y]=interseccion(f1,f2,g1,g2,x1,x2)
  %pendientes
  mf=(f2-f1)/(x2-x1);
  mg=(g2-g1)/(x2-x1);
  %ordenada al origen
  bf=f1-mf*x1;
  bg=g1-mg*x1;
  %Interseccion
  x=(bf-bg)/(mg-mf);
  y=mf*x+bf;
end
