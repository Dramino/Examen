%% C�digo para la obtenci�n de la respuesta en el dominio del tiempo
% Programa para obtener la respuesta de un puente con el m�todo de Strommen

%% Preparativos y captura de datos
clear all
close all
% Captura de datos
datosPuente  = dlmread('Puente.txt','',1,0);
datosVientoR = dlmread('VientoRafaga.txt','',1,0);
% LLenado de datos
puente = Puente;
viento = Viento;

% Asigna los datos a cada objeto
puente.L    = datosPuente (:,2);
puente.Lexp = datosVientoR(:,1);
puente.m    = datosPuente (:,3);
puente.fiYP = datosPuente (:,4);
puente.fiZP = datosPuente (:,5);
puente.fiXP = datosPuente (:,6);
puente.fiY  = datosVientoR(:,4);
puente.fiZ  = datosVientoR(:,5);
puente.fiX  = datosVientoR(:,6);
puente.w    = [datosPuente(1,8); datosPuente(2,8)];
puente.c    = [datosPuente(1,7); datosPuente(2,7)];
% Punto a revisar
yr = datosVientoR(1,15);
zr = datosVientoR(1,16);
B  = datosVientoR(1,17);
% Geometr�a
puente.B  = datosVientoR(:,2);
puente.D  = datosVientoR(:,3);
% Coeficiente aoerdinm�micos
puente.cD  = datosVientoR(:,7);
puente.cL  = datosVientoR(:,8);
puente.dCd = datosVientoR(:,9);
puente.cdl = datosVientoR(:,10);
% Obtenci�n de masa modal
puente.mModal();
% Integral de fi
puente.integrarFi()
% Datos del viento
viento.rho=1.25;
viento.v=datosVientoR(1,11);
viento.h=datosVientoR(1,12);

% Turbulencia
puente.Iu=datosVientoR(1,13);
puente.Iv=datosVientoR(1,14);
puente.integrarFi();


%Amortiguamiento aerodin�mico
%Definici�n de integrales
intY   = 0;
intFiY = 0;
intZ   = 0;
intFiZ = 0;

%Buffeting

%% Integrales
% $$\int_{L_{exp}}{\phi^{2}_{y}\cdot D \cdot\bar{C_D}dx}$$
%%
% $$\int_{L_{exp}}{\phi^{2}_{z}\cdot B \cdot C^{'}_{L} \cdot D\cdot\bar{C_D}dx}$$
%%
% $$\int_{L}{\phi^{2}_{n}dx}$$
intY   = integrarS(puente.fiY.^2.*puente.D.* puente.cD,puente.Lexp);
intZ   = integrarS(puente.fiZ.^2.*(puente.B.*puente.cdl + puente.D.*puente.cD),puente.Lexp);
%Puente
intFiY = integrarS(puente.fiYP.^2.,puente.L);
intFiZ = integrarS(puente.fiZP.^2.,puente.L);

%% Amortiguamiento aerodin�mico
% $${\zeta }_{ae_y}=-\frac{\rho}{2\tilde{m}_y}\cdot\bar{C}_{D}D\left( \frac{V}{\omega}\right)\frac{\int_{L_{exp}}{\phi_{y}^{2}dx}}{\int_L{\phi_{y}^{2}dx}}$$
%%
% $${\zeta }_{ae_y}=-\frac{\rho B^2}{4\tilde{m}_y}\cdot\left(C^{'}_{L}+ \bar{C}_{D}\frac{D}{B}\right)\left( \frac{V}{B\omega}\right)\frac{\int_{L_{exp}}{\phi_{y}^{2}dx}}{\int_L{\phi_{y}^{2}dx}}$$

zetaAeY = -(viento.rho*viento.v) /(2* puente.w(1)* puente.m1) * intY / intFiY ;
zetaAeZ = -(viento.rho*viento.v) /(4* puente.w(2)* puente.m2) * intZ / intFiZ ;
zetaAeM = 0

%% Funci�n de respuesta mec�nica H
% $${{\hat{H}}_{n}}\left( \omega  \right)={{\left[ 1-{{\left( \frac{\omega }{{{\omega }_{n}}} \right)}^{2}}+2i\left( {{\zeta }_{y}}-{{\zeta }_{a{{e}_{n}}}} \right)\cdot \frac{\omega }{{{\omega }_{n}}} \right]}^{-1}}$$
Hy = abs (1 - ( viento.w ./ puente.w(1) ) .^2 + 2*1i *(puente.c(1) - zetaAeY) *(( viento.w ./ puente.w(1) ))) .^( -1) ;
Hz = abs (1 - ( viento.w ./ puente.w(2) ) .^2 + 2*1i *(puente.c(2) - zetaAeZ) *(( viento.w ./ puente.w(2) ))) .^( -1) ;

%% Funcion de densidad espectral de Kaimal
%  Se obtiene apartir de la funci�n dentro de los m�todos de la clase Viento
[sNu,sNw] = viento.densidadEspectralKaimal();

%% Funci�n de acpentancia conjunta normalizada
% $$\hat{C}o_{nn}=\exp\left(-c_ny\cdot\frac{\omega\Delta S}{2\pi V}\right)$$

%%
% $$J_{y}^{2} = \int\int\phi_y(x_1)\cdot \phi_y(x_1)\cdot  \left\{ \left( 2\frac{D}{B}\bar{C}_D I_u\right)^2S_{uu}\hat{C}o_{uu}+ \left[\left(\frac{D}{B}C^{'}{L}-\bar{C}_L \right)I_w\right]^2 S_{ww}\hat{C}o_{ww}\right\} \,dx_1\,dx_2 $$
%%
% $$J_{z}^{2} = \int\int\phi_z(x_1)\cdot \phi_z(x_1)\cdot  \left\{ \left( \bar{C}_L I_u\right)^2 S_{uu}\hat{C}o_{uu}+ \left[\left(C^{'}_{L}+\frac{D}{B}\bar{C}_{D} \right)I_w\right]^2 S_{ww}\hat{C}o_{ww}\right\} \,dx_1\,dx_2 $$
%Co espectro - Constantes
coU = zeros(1 ,length(viento.w));
coW = zeros(1 ,length(viento.w));
cuY = 9;
cwY = 6;

%Funci�n de aceptancia constantes
jY = 0;
jZ = 0;
jm = 0;
%Expansi�n y prolongaci�n de datos
div = 4;
puente.D     = puente.prolong ( puente.D   , div );
puente.B     = puente.prolong ( puente.B   , div );
puente.cD    = puente.prolong ( puente.cD  , div );
puente.cL    = puente.prolong ( puente.cL  , div );
puente.dCd   = puente.prolong ( puente.dCd  , div );
puente.cdl   = puente.prolong ( puente.cdl  , div );
puente.fiY   = puente.prolong ( puente.fiY , div );
puente.fiZ   = puente.prolong ( puente.fiZ , div );
puente.Lexp  = puente.expand  ( puente.Lexp, div );
puente.L     = puente.expand  ( puente.L   , div );


n = length(puente.Lexp);
x=cumsum(puente.Lexp);
jY=0;
jZ=0;
% funci�n de acepantcia conjunta
for i=1:n
  for j=1:n
        Dx=abs(x(i)-x(j));
        coU = exp (- cuY* ( viento.w*Dx ) /(2 * pi * viento.v));
        coW = exp (- cwY* ( viento.w*Dx ) /(2 * pi * viento.v));
    jY  = jY + puente.fiY(i) * puente.fiY(j) .* ( (((2* puente.Iu )/ B) .^2* puente.D(i) * puente.cD(i) * puente.D(j) * puente.cD(j) ) .* coU.* sNu + ( (puente.D(i)/puente.B(i)*puente.dCd(i)-puente.cL(i)) * puente.Iv ) .^2 .* coW.* sNw ) .* puente.L(i) .* puente.L(j);
        jZ  = jZ + puente.fiZ(i) * puente.fiZ(j) .* ( (2* puente.cL(i)* puente.Iu) ^2 .* coU.* sNu + (( puente.Iv^2*( puente.cdl(i) +( puente.D(i) * puente.cD(i) )./B) *(puente.cdl(j) +( puente.D(j) * puente.cD(j) )./ B)) .*coW.* sNw )) * puente.L(i) * puente.L(j) ;
    %jm  = jm + puente.fim(i) * puente.fim(j) .* ( (2* puente.cM(i)* puente.Iu) ^2 .* coU.* sNu + (( puente.Iv^2*( puente.cdm(i)*(puente.cdm(j))))) .*coW.* sNw )) * puente.L(i) * puente.L(j) ;
  end
end

% Normalizaci�n de la funci�n de aceptancia
jNormY = jY /( [puente.intFi1] ) ^2;
jNormZ = jZ /( [puente.intFi2] ) ^2;
%jNormM = jM /( [puente.intFi3] ) ^2;

%% Espectro de respuesta
% $${{S}_{{{r}_{n}}}}\left( {{x}_{r}},\omega  \right)={{\left[ \frac{\rho {{V}^{2}}B}{2}\cdot {{{\hat{J}}}_{n}}\left( \omega  \right) \right]}^{2}}\hat{H}_{n}^{y}$$
espResDivFiY = (( viento.rho *viento.v ^2.* B) /(2* puente.m1 * puente.w(1)^2) ) ^2 * Hy.^2.* jNormY;
espResDivFiZ = (( viento.rho *viento.v ^2.* B) /(2* puente.m2 * puente.w(2)^2) ) ^2 * Hz.^2.* jNormZ;
%Respuesta al extremo del puente
espResY = espResDivFiY*yr^2;
espResZ = espResDivFiZ*zr^2;




%% Desviaci�n est�ndar
% $$\sigma_n(x_r)=\left|\phi_n(x_r)\right|\cdot\left(\frac{\rho B^3}{2\tilde{m}_n}\right)^2\cdot\left[\int_{0}^{\infty}\left|\hat{H}_n(\omega)\right|^2\cdot\hat{J}_n^2(\omega)d\omega\right]^{1/2}$$
Nomega = length ( viento.w);
Domega = ( viento.w ( Nomega ) - viento.w(1) ) /( Nomega -1) ;
intSy  = integrar(Hy.^2.*jNormY,Domega);
intSz  = integrar(Hz.^2.*jNormZ,Domega);

sRy = abs(yr)*viento.rho *viento.v^2*B/(2*puente.m1*puente.w(1)^2)*sqrt(intSy);
sRz = abs(zr)*viento.rho *viento.v^2*B/(2*puente.m2*puente.w(2)^2)*sqrt(intSz);
%% Respuesta en el dominio del tiempo
% $$x(t)=\sum\limits_{k=1}^{N}{{{c}_{k}}\cos \left( {{\omega}_{k}}t+{{\theta }_{k}} \right)}$$
%Divisions de la frecuencia del viento
Nw = length ( viento.w );
dw = ( viento.w( Nw ) - viento.w(1)) / Nw ;
%Simulaci�n de 10 minutos
t   = linspace (0 ,600 ,600) ;
ry  = zeros ( 1 , length(t) );
rz  = zeros ( 1 , length(t) );
for i = 1: Nw
  ry = ry+ sqrt(2* espResY(i)* dw ) .* cos ( viento.w(i)*t - rand (1) *2* pi);
  rz = rz+ sqrt(2* espResZ(i)* dw ) .* cos ( viento.w(i)*t - rand (1) *2* pi);
end

%Factor pico
maxY = max(abs(ry));
maxZ = max(abs(rz));
kpY  = maxY / sRy;
kpZ  = maxZ / sRz;

%Gr�fica de la respuesta
subplot (2 ,1 ,1)
plot (t , ry,'b')
grid
title('Desplazamientos horizontales')
xlabel ( '$T$ [s] ','Interpreter','latex')
ylabel ( '$r_y$ [m] ','Interpreter','latex')

hold all

subplot (2 ,1 ,2)
plot (t , rz,'r' )
grid
xlabel ( '$T$ [s] ','Interpreter','latex')
ylabel ( '$r_z$ [m] ','Interpreter','latex')
title('Desplazamientos verticales')

%% Resultados
archivo= fopen ('Resultado.txt', 'w');
fprintf(archivo,'%s\r\n',datetime);
fprintf(archivo,'Desplazamiento máximo horizontal: %6.4f m\r\n',maxY);
fprintf(archivo,'Desplazamiento máximo vertical: %6.4f m\r\n',maxZ);
fprintf(archivo,'Desviación estándar horizontal: %6.4f m\r\n',sRy);
fprintf(archivo,'Desviación estándar vertical: %6.4f m\r\n',sRz);
fprintf(archivo,'Factor pico horizontal: %6.4f\r\n',kpY);
fprintf(archivo,'Factor pico vertical: %6.4f\r\n\r\n',kpZ);
fprintf(archivo,'%12s %12s %12s\r\n','T[s]', 'Desp Hor [m]', 'Desp Ver [m]');
fprintf(archivo,'%12.4f %12.4f %12.4f\r\n',[t;ry;rz]);
fclose(archivo);
