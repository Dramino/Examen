%% Integral por el m�todo de simpson para distancia equidistantes
% $$\int_{a}^{b}f(x)dx\approx\frac{h}{3}\left[ f(x_0)+2\sum_{i=1}^{(n/2-1)}f(x_{2j})+4\sum_{i=1}^{n/2}f(x_{2j-1}+f(x_n)\right]$$
function int=integrar(fx,Dx)
  %integral de una serie de datos por medio del método de simpson
  n=length(fx);
    int=0;
    for i=1 : n
    if i==1 || i==n
      int=int+fx(i);
    else
      if mod(i,2)==0
        int=int+4*fx(i);
      else
        int=int+2*fx(i);
      end
    end
    end
    int=int*Dx/3;
end
