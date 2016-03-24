function u = Tridiagonal( a, b, c, r )

n = length(r);

u = zeros(n,1);   
gam = zeros(n,1);   

bet = b(1);
u(1) = r(1)/bet;

for j=2:n
    gam(j) = c(j-1)/bet;
    bet = b(j) - a(j) * gam(j);
    u(j) = (r(j)-a(j)*u(j-1))/bet;
end

for i=1:n-1
    j = n - i;
    u(j) = u(j) - gam(j+1) * u(j+1);
end
    