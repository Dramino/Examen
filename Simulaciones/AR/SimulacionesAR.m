% SimulacionesAR.m
%% Simulaciones por el método AR
%% Importar datos
clear
close
proVie=3; proSim=4; proLog=4; propiedades=proVie+proSim+proLog;

dato= importar(propiedades);

V        = dato(1);
z        = dato(2);
z0       = dato(3);
t        = dato(4);
ti       = dato(5);
ds       = dato(6);
N        = dato(7);
cn       = dato(8);
graficar = dato(9);
exportar = dato(10);
nodo     = dato(11);

%% Inicio del programa
dx = 0;
U  = 0.4*V/(log(z/z0));
P  = 4;
x  = cumsum(zeros(N,1)+ds);
n  = floor(t/ti);

% Revivisón de coherencia
if cn
  cnu=9; cnw=6;
else
  cnu=0; cnw=0;
end

Au=zeros(P*N);
Aw=zeros(P*N);
Ru=zeros(N);
Rw=zeros(N);
for p=1:P
  R=zeros(N);
  for i=1:N
    for j=i:N
      dx     = abs(x(i)-x(j));
      sNu    = @(f)U^2*105*(z/V)./(1+33*f/(2*pi)*z/V).^(5/3).*exp(-cnu*dx*f/(2*pi)/V).*cos(2*pi*f/(2*pi)*(p-1)*ti);
      sNw    = @(f)U^2*2*(z/V)./(1+5.3*f/(2*pi)*z/V).^(5/3).*exp(-cnw*dx*f/(2*pi)/V).*cos(2*pi*f/(2*pi)*(p-1)*ti);

      Ru(i,j) = integral(sNu,0,50);
      Ru(j,i) = Ru(i,j);

      Rw(i,j) = integral(sNw,0,50);
      Rw(j,i) = Rw(i,j);
    end
    end
    for k=1:P-p+1
        Au(((k-1)*N+1):(k*N),((k+p-1-1)*N+1):((k+p-1)*N))=Ru;
    Aw(((k-1)*N+1):(k*N),((k+p-1-1)*N+1):((k+p-1)*N))=Rw;
    end
end
% Matriz simétrica
for i=1:P*N
    for j=1:i
        Au(i,j)=Au(j,i);
    Aw(i,j)=Aw(j,i);
    end
end

% Correlación sin P
Ru=zeros(N);
Rw=zeros(N);
for i=1:N
  for j=i:N
    dx     = abs(x(i)-x(j));
    sNu    = @(f)U^2*105*(z/V)./(1+33*f/(2*pi)*z/V).^(5/3).*exp(-cnu*dx/(2*pi)*f/V).*cos(2*pi*f/(2*pi)*P*ti);
    sNw    = @(f)U^2*2*(z/V)./(1+5.3*f/(2*pi)*z/V).^(5/3).*exp(-cnw*dx/(2*pi)*f/V).*cos(2*pi*f/(2*pi)*P*ti);

    Ru(i,j) = integral(sNu,0,50);
    Ru(j,i) = Ru(i,j);

    Rw(i,j) = integral(sNw,0,50);
    Rw(j,i) = Rw(i,j);
  end
end
Bu  = Au(1:N,(N+1):(N*P));  Bw  = Aw(1:N,(N+1):(N*P));
Bu  = [Bu,Ru];              Bw  = [Bw,Rw];
Bu  = Bu';                  Bw  = Bw';
Xu  = Au\Bu;                Xw  = Aw\Bw;
R0u = Au(1:N,1:N);          R0w = Aw(1:N,1:N);
RNu = R0u;                  RNw = R0w;

for i=1:P
  RNu = RNu-(Xu(((i-1)*N+1):(i*N),1:N))'*Bu(((i-1)*N+1):(i*N),1:N);
  RNw = RNw-(Xw(((i-1)*N+1):(i*N),1:N))'*Bw(((i-1)*N+1):(i*N),1:N);
end

Lu=chol(RNu);    Lw=chol(RNw);
Lu=Lu';          Lw=Lw';
au=zeros(N,n);   aw=zeros(N,n);
for i=1:N
    au(i,:)=normrnd(0,1,1,n);
  aw(i,:)=normrnd(0,1,1,n);
end

Vu(1:N,1) = Lu*au(:,1);
Vw(1:N,1) = Lw*aw(:,1);

for i=2:p
    Vu(1:N,i) = Lu*au(:,i); Vw(1:N,i) = Lw*aw(:,i);
    for j=1:i-1
        Vu(1:N,i) = (Xu(((j-1)*N+1):(j*N),:))'*Vu(1:N,(i-j))+Vu(1:N,i);
    Vw(1:N,i) = (Xw(((j-1)*N+1):(j*N),:))'*Vw(1:N,(i-j))+Vw(1:N,i);
    end
end

for i=(P+1):n
    Vu(1:N,i) = Lu*au(:,i);  Vw(1:N,i) = Lw*aw(:,i);
    for j=1:P
        Vu(1:N,i)=(Xu(((j-1)*N+1):(j*N),:))'*Vu(1:N,(i-j))+Vu(1:N,i);
    Vw(1:N,i)=(Xw(((j-1)*N+1):(j*N),:))'*Vw(1:N,(i-j))+Vw(1:N,i);
    end
end

t=(1:n)*ti;

if graficar; crearGrafica(t,Vu(nodo,:),Vw(nodo,:));end
if exportar; crearTxt(Vu(nodo,:),Vw(nodo,:));end



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
