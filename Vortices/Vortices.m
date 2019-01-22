%Vortices.m
%% Respuesta del puente ante efecto de vórtices
% Obtiene la respuesta en el dominio del tiempo
clear all
close all

%Captura de datos
datosPuente  = dlmread('Puente.txt','',1,0);
datosVientoV = dlmread('VientoVortices.txt','',1,0);
%LLenado de datos
puente = Puente;
viento = Viento;

%Asigna los datos a cada objeto
puente.L    = datosPuente (:,2);
puente.Lexp = datosVientoV(:,1);
puente.m    = datosPuente (:,3);
puente.fiYP = datosPuente (:,4);
puente.fiZP = datosPuente (:,5);
puente.fiXP = datosPuente (:,6);
puente.fiZ  = datosVientoV(:,2);
puente.fiX  = datosVientoV(:,3);
puente.w    = [datosPuente(1,8); datosPuente(2,8)];
puente.c    = [datosPuente(1,7); datosPuente(2,7)];
puente.integrarFi();
puente.mModal();

%Propiedades del viento
viento.Kaz = datosVientoV(1,4);
viento.az  = datosVientoV(1,5);

%Datos de espectro de carga
viento.lambdaZ = datosVientoV(1,6);
viento.sQz     = datosVientoV(1,7);
viento.bZ      = datosVientoV(1,8);

%Datos del viento
viento.rho = 1.25;
viento.v   = datosVientoV(1,9);
puente.St  = datosVientoV(1,10);

%Valor del modo al extremo del puente en direccion Z
fiI = datosVientoV(1,11);
wZ  = puente.w(2);
mZ  = puente.m2;
D   = datosVientoV(1,13);
B   = datosVientoV(1,12);
%% Desviación estándar
% $$\sigma_{r_z}=\hat{\sigma}_{r_z}\cdot{a_z}D$$
%%
% $${{\hat{\sigma }}_{{{r}_{\theta }}}}={{\left\{ \frac{1-\hat{\zeta }}{2}+{{\left[ {{\left( \frac{1-\hat{\zeta }}{2} \right)}^{2}}+{{{\hat{\beta }}}^{2}} \right]}^{1/2}} \right\}}^{1/2}}$$
%% 
% $$\hat{\zeta }=\frac{4{{{\tilde{m}}}_{\theta }}}{\rho {{B}^{4}}}\cdot \frac{{{\zeta }_{\theta }}}{{{K}_{{{a}_{\theta }}}}}\cdot \frac{\int\limits_{L}{\phi _{\theta }^{2}dx}}{\int\limits_{L\exp }{\phi _{\theta }^{2}dx}}$$
%%
% $$\hat{\beta }=\frac{\left| {{\phi }_{\theta }}\left( {{x}_{r}} \right)\right|}{{{2}^{5/2}}{{\pi }^{7/2}}}\cdot {{\left( \frac{\rho{{D}^{5}}}{{{{\tilde{m}}}_{\theta }}\int\limits_{L}{\phi _{\theta }^{2}dx}}\cdot \frac{\lambda }{{{b}_{\theta }}{{K}_{{{a}_{\theta }}}}} \right)}^{1/2}}\cdot \frac{{{{\hat{\sigma }}}_{{{q}_{\theta }}}}}{S{{t}^{2}}}\cdot \frac{{{g}_{\theta }}}{{{a}_{\theta }}}$$
puente.intFi2=integrarS(puente.fiXP.^2+puente.fiZP.^2,puente.L);
intVor = integrarS(( puente.fiX.^2 + puente.fiZ.^2),puente.Lexp);
Vr     = D* wZ /(2* pi* puente.St ); %Velocidad resonante
V      = D*puente.w(2)/(2*pi*puente.St);
z      = ((4* puente.m2 * puente.c(2) ) /( viento.rho *B^2* viento.Kaz )) *( intVor / puente.intFi2 );
b      = fiI /(2^(5/2) * pi ^(7/4) ) * (( viento.rho *D^3* viento.lambdaZ ) /( puente.m2 * viento.bZ * viento.Kaz * puente.intFi2) )^(1/2) * ( viento.sigmaQz /( puente.St ^2* viento.az ) ) * (V / Vr ) ^(3/2) * exp ( -(1/2) *((1 -( Vr /V)) / viento.bZ )^2) ;
sAux   = ((1 - z) /2 + (((1 - z ) /2) ^2 + b ^2) ^(1/2) ) ^(1/2) ;
sRz    = sAux * viento.az*D;

%% Coeficiente de amortiguamiento aerodinámico
% $${{\zeta }_{a{{e}_{z}}}}=\frac{\rho {{B}^{2}}}{4{{m}_{z}}}\cdot{{K}_{{{a}_{z}}}}\left[ 1-{{\left( \frac{{{\sigma }_{z}}}{{{a}_{z}}D} \right)}^{2}} \right]\cdot \frac{\int\limits_{L\exp }{\phi _{z}^{2}dx}}{\int\limits_{L}{\phi _{z}^{2}dx}}$$
zAe=((viento.rho*B^2*viento.Kaz)/(4*mZ)).*((1-(sRz./(viento.az*D)).^2)*(puente.intFi2/intVor));

%% Función de respuesta H de frecuencias
% $${{\hat{H}}_{n}}\left( \omega  \right)={{\left[ 1-{{\left( \frac{\omega }{{{\omega }_{n}}} \right)}^{2}}+2i\left( {{\zeta }_{y}}-{{\zeta }_{a{{e}_{n}}}} \right)\cdot \frac{\omega }{{{\omega }_{n}}} \right]}^{-1}}$$
H = abs (1 - ( viento.w ./ wZ ) .^2 + 2*1i *( puente.c(2) - zAe ) *(( viento.w./ wZ ))).^( -1) ;

%% Espectro de carga
% $${{S}_{\hat{Q}n}}(\omega )=2\lambda D\cdot \frac{{{\left(
% \frac{1}{2}\rho {{V}^{2}}B{{{\hat{\sigma }}}_{{{q}_{z}}}}
% \right)}^{2}}}{\sqrt{\pi }{{b}_{s}}{{\omega }_{s}}}\cdot {{e}^{-{{\left( \frac{1-\omega /{{\omega }_{s}}}{{{b}_{z}}} \right)}^{2}}}}\frac{\int\limits_{L\exp }{\phi _{n}^{2}\left( x \right)dx}}{{{\left( \omega _{n}^{2}{{{\tilde{M}}}_{n}} \right)}^{2}}}$$

wS = (2* pi* puente.St *Vr)/D ; %Frecuencia resonante
sQ = 2* viento.lambdaZ * D * (((0.5* viento.rho *Vr ^2* B* viento.sQz) ^2) /( sqrt (pi)* wS * viento.bZ )) * exp( -(((1 - viento.w ./ wS ) ./ viento.bZ ) .^2)) .* intVor ;

%% Espectro de respuesta
% $${{S}_{{{r}_{n}}}}({{x}_{r}},\omega )=\phi _{i}^{2}\left( {{x}_{r}}\right){{\left| \hat{H}\left( \omega  \right) \right|}^{2}}\cdot{{S}_{{{{\tilde{Q}}}_{n}}}}\left( \omega  \right)$$
espRes = (( fiI^2* H .^2) ./(( wZ ^2* mZ * puente.intFi2) .^2) ) .* sQ ;

%% Respuesta al extremo del puente en el dominio del tiempo
% $$x(t)=\sum\limits_{k=1}^{N}{{{c}_{k}}\cos \left( {{\omega}_{k}}t+{{\theta }_{k}} \right)}$$
t  = linspace (0 ,600 ,600) ;
n  = length([viento.w]);
dw = ( viento.w (n)- viento.w(1)) / n ;
rz = zeros (1 ,length(t));
for i = 1: n
	rz = rz + sqrt (2* espRes(i)* dw ) .* cos ( viento.w(i)*t - rand (1) *2* pi);
end

%Valor máximo
maxZ = max (abs (rz));
kpz  = maxZ / sRz;

%Gráfica de la respuesta en el tiempo
plot (t , rz,'b')
grid
xlabel ( '$T [s]$','Interpreter','latex','LineWidth',1.25)
ylabel ( '$r_z [m]$','Interpreter','latex','LineWidth',1.25)
title('Desplazamientos verticales')

%Impresión de resultados
archivo= fopen ('Resultado.txt', 'w');
fprintf(archivo,'%s\r\n',datetime);
fprintf(archivo,'Desplazamiento máximo vertical: %6.4f m\r\n',maxZ);
fprintf(archivo,'Desviación estándar vertical: %6.4f m\r\n',sRz);
fprintf(archivo,'Factor pico vertical: %6.4f\r\n\r\n',kpz);
fprintf(archivo,'%12s %12s \r\n','T[s]','Desp Ver [m]');
fprintf(archivo,'%12.4f %12.4f \r\n',[t;rz]);
fclose(archivo);
