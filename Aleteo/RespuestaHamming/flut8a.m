function dx=flut8a(t,x,C,B,U,V,dt)


if not(mod(t,dt)==0)
	n0 = floor(t/dt)+1;
    if n0==length(U); n0=length(U)-1; end
	n1 = n0+1;
	t0 = dt*(n0-1);
	t1 = t0+dt;
	u  = [U(n0) U(n1)];
	v  = [V(n0) V(n1)];
	T  = [t0 t1];
    F1 = griddedInterpolant(T(:),u(:),'linear');
    F2 = griddedInterpolant(T(:),v(:),'linear');
	U  = F1(t);
	V  = F2(t);
else
	n  = floor(t/dt)+1;
    if n > length(U) ; n=length(U); end
	U  = U(n);
	V  = V(n);
end
dx(1)=x(3);
dx(2)=x(4);
dx(3)=C(1)*x(1)+C(2)*x(2)+C(3)*x(3)+C(4)*x(4)+B(1)*U + B(2)*V;
dx(4)=C(5)*x(1)+C(6)*x(2)+C(7)*x(3)+C(8)*x(4)+B(3)*U + B(4)*V;
