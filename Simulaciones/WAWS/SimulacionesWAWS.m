%% Programa para la obtención de simulaciones coherente e incoherente con el método WAWS

%% Importación de datos
clear
close
proVie = 3;
proSim = 4;
proLog = 4;
propiedades = proVie + proSim + proLog;
datos       = importar(propiedades);


v            = datos(1);
h            = datos(2);
z0           = datos(3);
t            = datos(4);
dt           = datos(5);
separacion   = datos(6);
simulaciones = datos(7);
cn           = datos(8);
graficar     = datos(9);
exportar     = datos(10);
nodo         = datos(11);

%Revisa si hay coherencia
if cn==true
  cu=9;  cw=6;
else
  cu=0;  cw=0;
end

[sNu,sNw,w] = densidadKaimal(h,v,z0);
intTiempo = [0 t];
N         = t/dt;
dw        = (w(length(w)) - w(1)) / length(w);
u         = WindSimMult(sNu,simulaciones,separacion,v,dw,intTiempo,N,cu);
w         = WindSimMult(sNw,simulaciones,separacion,v,dw,intTiempo,N,cw);
T         = [0, linspace(dt,t,N)];
if graficar; crearGrafica(T,u(nodo,:),w(nodo,:));end
if exportar; crearTxt(u(nodo,:),w(nodo,:));end

%% Función de densidad de kaimal
function [sNu,sNw,w]=densidadKaimal(z,V,z0)
      %Obtiene la funciÃ³n de densidad espectral en todo el dominio w

  U   =0.4*V/(log(z/z0));
  %Definicion de variables
  w   = linspace(0,10,2000);
  sNu =U^2*105*(z/V)./(1+33*w/(2*pi)*z/V).^(5/3);
  sNw =U^2*2*(z/V)./(1+5.3*w/(2*pi)*z/V).^(5/3);
end

%% Función para las simulaciones
function [x]=WindSimMult(Sw,numSim,separacion,velocidad,Dw,tiempo,Ndatos,coherencia)

  numFreq  = 2000;
  dx       = separacion;
  t0       = tiempo(1);
  tf       = tiempo(2);
  t        = [t0,linspace(dx,tf,Ndatos)];
  w0       = 0;
  w        = w0;
  x        = zeros(numSim,length(t));
  fi       = 2 * pi * rand(numSim,numFreq);


  if not(coherencia==0)
    % Con coherencia
    for ii = 1:numFreq
      cx = coherencia;
      Co = exp( -( cx *w* dx) / (2 * pi * velocidad));
      %Cholesky decomposition
      G = choleskySimplify(numSim,Co);

      a = sqrt(2 * Sw(ii) * Dw);

      for jj=1:numSim
        for k=1:jj
          x(jj,:) = x(jj,:) + G(jj,k) .* a .* cos( w .* t + fi(k,ii));
        end
      end
      w = w + Dw;

        end
  else
    %Sin coherencia
        for ii=1:numSim

      for jj = 1: numFreq
        x(ii,:) = x(ii,:) + sqrt(2* Sw(jj)* Dw ) .* cos ( w.*t + rand (1) *2* pi);
        w       = w + Dw;
      end
    end
  end
end


function [G]=choleskySimplify(matDimM,Coh)
  C = Coh;
  G = zeros(matDimM);
  for j=1:matDimM
    for m = 1:j
      if m == 1 && m<=j
        G(j,m) = C ^ (abs(j-m));
      else
        G(j,m) = C ^ (abs(j-m)) * sqrt(1-C ^ 2);
      end
    end
  end
end




%% Función para importar datos
function mDatos=importar(numDatos)

  [archivo,direccion] = uigetfile({'*.txt'},'File Selector');

  dirCom     = strcat(direccion,archivo);
  archivoID  = fopen(dirCom,'r');
  datos      = textscan(archivoID, '%s','Delimiter', '\n');
  matriz     = datos{1};

  fclose(archivoID);

  varConst  = zeros(1,numDatos);
  revConst  = false;
  numConst  = 0;
  for i = 1 : length(matriz)
    if revConst == false
      constante = strsplit(matriz{i},'=');
      if length (constante)== 2
        num                = str2double(strsplit(constante{2}));
        numConst           = numConst + 1;
        varConst(numConst) = num(2);
        if numConst == numDatos, revConst=true; end
      end

    end
  end
  mDatos=varConst;
end

%% Función para exportar
function crearTxt(u,w)
  fileID = fopen('vel.txt','w');
  fprintf(fileID,'%10s %10s\r\n','u','w');
  fprintf(fileID,'%10.4f %10.4f\r\n',[u(:) w(:)]);
  fclose(fileID);
end

%% Función para graficar
function crearGrafica(T,u,w)
  subplot(2,1,1)
  plot(T,u,'b')
  xlabel('T [s]')
  ylabel('V [m/s]')
  title('Velocidad fluctuante u');

  subplot(2,1,2)
  plot(T,w,'b')
  xlabel('T [s]')
  ylabel('V [m/s]')
  title('Velocidad fluctuante w');
end
