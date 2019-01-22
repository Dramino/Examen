%% Función para obtener los las derivadas aerodinámicas
%derAer.m
function [H1,H2,H3,H4,A1,A2,A3,A4]=derAer(vMin,vMax,Datos)
  if vMin<=0, vMin=0.01; end
  v=linspace(vMin,vMax,Datos);

  x  = 0.5./v;
  j0 = besselj(0,x);
  j1 = besselj(1,x);
  y0 = bessely(0,x);
  y1 = bessely(1,x);

  G = -(j1.*j0+y1.*y0)./((j1+y0).^2+(y1-j0).^2);
  F = (j1.*(j1+y0)+y1.*(y1-j0))./((j1+y0).^2+(y1-j0).^2);

  H1 = -2*pi*v.*F;
  H2 = pi/2*(1+F+4*G.*v).*v;
  H3 = 2*pi*(F.*v-G/4).*v;
  H4 = pi/2*(1+4*G.*v);

  A1 = -pi/2*F.*v;
  A2 = -pi/8*(1-F-4*G.*v).*v;
  A3 = pi/2*(F.*v-G/4).*v;
  A4 = pi/2*G.*v;
end
