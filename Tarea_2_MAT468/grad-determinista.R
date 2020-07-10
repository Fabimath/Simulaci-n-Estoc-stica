library(numDeriv)                           

# Mostramos un mapa del dominio de h(x,y)
h = function(x, y)  (x*sin(20*y) + y*sin(20*x))^2*cosh(sin(10*x)*x)+ (x*cos(10*y) - y*sin(10*x))^2*cosh(cos(20*y)*y)
a = 1
x = seq(from=-a, to=a, length.out=1000)
y = seq(from=-a, to=a, length.out=1000)
z = outer(x, y, h)

# Tomamos nuestra estimación inicial
theta_0 = c(0.4, 0.5)
alfa = 0.1
# Número de iteraciones
k = 150

f = function(x)  (x[1]*sin(20*x[2]) + x[2]*sin(20*x[1]))^2*cosh(sin(10*x[1])*x[1])+ (x[1]*cos(10*x[2]) - x[2]*sin(10*x[1]))^2*cosh(cos(20*x[2])*x[2])


thetas = matrix(NA, ncol=2, nrow=k+1)      
colnames(thetas) = c('x', 'y')
thetas[1, ] = theta_0                      

for (i in 2:(k+1))
  thetas[i, ] = thetas[i-1, ] - alfa * grad(func=f, x=thetas[i-1, ])

image(x=x, y=y, z=z)
contour(x=x, y=y, z=z, add=TRUE, nlevels=40, col=gray(0.5))
points(theta_0[1],theta_0[2],col='blue',lwd = 3)
points(thetas[,1],thetas[,2],col='green',type='l')
points(thetas[k,1],thetas[k,2],col='red',lwd = 3)
