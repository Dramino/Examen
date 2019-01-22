classdef Puente<handle
   
    
    properties
        %Formas modales expuestos
        fiX = [];
        fiY = [];
        fiZ = [];
		%Formas modales totales del puente
		fiXP = [];
        fiYP = [];
        fiZP = [];
        %Longitudes caracteristicas del puente
        %B=ancho D=alto L=Longitud de cada elemento m=masa de cada elemento
        %Lexp=Longitud expuesta
        B    = [];
        D    = [];
        L    = [];
		Lexp = [];
        %Propiedade dinámicas del viento
        m    = [];  %masa de cada elemento
        m1   = 0 ;  %masa modal efectiva del modo1
        m2   = 0 ;  %masa modal efectiva del modo2
        w    = [];  %Frecuencias de cada modo
        c    = [];  %Coeficiente de amortiguamiento de cada modo
        cD   = [];  %Coeficiente de arrastre
        cL   = [];  %coelficiente de levante
		dCd  = [];  %Derivada del coeficiente de arrastre
        cdl  = [];  %Derivada del coelficiente de levante
        St   = 0 ;  %Número de Strouhal
		Iu   = 0 ;  %Intensidad de turvulencia en direecion x
		Iv   = 0 ;  %Intensidad de turvulencia en direecion y
        
        %Integrales de las formas modales al cuadrado
        intFi1  = 0;  %Integral del modo 1 direccion y
        intFi2  = 0;  %integral del modo 2 direccion x z
    end
    
    methods
        function mModal(obj)
            %Obtiene la masa modal efectiva para las dos direcciones
            n       = length([obj.L]);
            intM1   = 0;
            intPhi1 = 0;
            intM2   = 0;
            intPhi2 = 0;
            %Datos del puente
            %integral sobre L
           
            intM1   = integrarS([obj.fiYP].^2.* [obj.m],[obj.L]) ;
            intPhi1 = integrarS([obj.fiYP].^2,[obj.L]);
            intPhi2 = integrarS([obj.fiXP].^2+[obj.fiZP].^2,[obj.L]);
            intM2   = integrarS(([obj.fiXP].^2+[obj.fiZP].^2).*[obj.m],[obj.L]);
            
			[obj.m1] = intM1/intPhi1;
            [obj.m2] = intM2/intPhi2;
        end
        
		function integrarFi(obj)
            n = length([obj.L]);
            for i = 1 : n
				% Longtiud total
                [obj.intFi1] = integrarS([obj.fiYP].^2,[obj.L]);
                [obj.intFi2] = integrarS([obj.fiXP].^2+[obj.fiZP].^2,[obj.L]);
            end
        end
 
        function [ VarNueva ] = prolong (obj,var , div )
        % aumenta la cantidad de datos
            n            = length (var);
            VarNueva     = zeros (0, (n -1) * div +1) ;
            VarNueva (1) = var (1) ;
            for i = 0:n -2
                delta = ( var (i +2) - var (i +1) )/ div ;
                for j = 1: div
                    VarNueva ( div *i +j +1) = VarNueva ( div *i+j ) + delta ;
                end
            end
        end
        
        function [ VarNueva ] = expand (obj, var , div )
            n                  = length (var);
            VarNueva           = zeros (0, (n - 1) * div + 1) ;
            VarNueva (2: div ) = (var (1) + var (2) * 0.5) / div ;
            VarNueva (1)       = 0.5 * VarNueva (2) ;
            for i = 1 : n -3
                delta                               = ( var (i +1) + var (i +2) ) /(2* div );
                VarNueva ( div *i +1: div *i+ div ) = delta ;
                VarNueva ( div *i + div )           = ( delta + ( var (i +2) + var (i +3) ) /(2* div )) /2;
            end
            VarNueva ((n -2) * div +1:( n -1) * div ) = ( var (n -1) *0.5+ var (n)) / div ;
            VarNueva ((n -1) * div +1)                = 0.5* VarNueva ((n -1) * div ) ;
            VarNueva ( div *(n -3) + div )            = ( delta + ( var (i +2) +2* var (i +3) ) /(2* div )) /2;
        end

    end
    
end

