classdef Viento

    properties
        rho = 1.25;  %Densidad del viento
        v   = 0;     %Velocidad media del viento
        h   = 0;     %Altura a revisar
        w   = linspace (0 ,20 ,2000) ;  %frecuencias del viento
        %Propiedades de desprendimiento de v�rtices
        Kaz     = 0.2;
        az      = 0.4;
        lambdaZ = 2;
        sQz     = 1;
        bZ      = 0.15;
        sigmaQz = 1
    end

    methods
        %% Funci�n de densidad de potencia espectral de Kaimal en u y w
        % $$\frac{S_n(\omega)}{\sigma_{m}^{2}}=\frac{A_n {^{x}L_n/V}}{(1+1.5\cdot A_n\omega {^{x}L_n/V)^{5/3}}}$$
        function [sNu,sNw]=densidadEspectralKaimal(obj)
            %Obtiene la funci�n de densidad espectral en todo el dominio w
            %para las dos direeciones
            %Constantes de la funcion
            Au  = 6.8/(2*pi);
            Aw  = 9.4/(2*pi);
            H   = [obj.h];
            xLu = 100*(H/10)^0.3;
            xLw = xLu/12;
            V   = [obj.v];
            %Definicion de variables
            w   = [obj.w];
            sNu = zeros(1,length(w));
            sNw = zeros(1,length(w));
            %Funci�n de densidad espectral
            sNu=Au.*xLu./(V *(1+1.5* Au.*w*(xLu/V )).^(5/3) );
            sNw=Aw.*xLw./(V *(1+1.5* Aw.*w*(xLw/V )).^(5/3) );
        end
    end

end
