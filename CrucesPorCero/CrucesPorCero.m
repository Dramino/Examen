% CrucesPorCero.m
% Calculo del cruces por cero y factor pico
dato = importar(2);
reg  = importarRegistro();

ti   = dato(1);
tf   = dato(2);
n    = length(reg);
t    = linspace(ti,tf,n);
tInt = zeros(n,1);
nu=0;
for i=1:n-1
    if reg(i)*reg(i+1)<0
      nu    = nu+1;
      m     = (reg(i+1)-reg(i))/(t(i+1)-t(i));
      b     = reg(i)-m*t(i);
      tInt(nu) = -b/m;
    end
end
grafTie   = tInt(1:nu);
grafCeros = zeros(length(grafTie));
kp        = sqrt(2*log(nu))+0.5772/(sqrt(2*log(nu)));

%Grafica
plot(t,reg,'b');
hold on
plot(grafTie,grafCeros,'or')

%Impresion de resultados
archivo= fopen ('Resultado.txt', 'w');
fprintf(archivo,'%s\r\n\r\n',datetime);
fprintf(archivo,'Total de veces que cruza por cero: %6.4f \r\n',nu);
fprintf(archivo,'Factor pico : %6.4f \r\n\r\n',kp);
fprintf(archivo,'Puntos en los que cruza por cero \r\n');
fprintf(archivo,'%12s \r\n','T [s]');
fprintf(archivo,'%12.4f \r\n',grafTie);
fclose(archivo);

function registro= importarRegistro()
  [archivo,direccion] = uigetfile({'*.txt'},'File Selector');
  dirCom   = strcat(direccion,archivo);
  registro = dlmread(dirCom,'',1,0);
end

%% Funcion de importar
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
