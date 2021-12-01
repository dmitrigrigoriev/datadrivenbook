library(pracma)
library(deSolve)
library(plotly)


a <- 1    # Thermal diffusivity constant
L <- 100  # Length of domain
N <- 1000 # Number of discretization points
dx <- L/N
x <- seq(-L/2, L/2, dx) # Define x domain

rhsHeat <- function(t, y, pars) {
  kappa <- pars[[1]]
  a <- pars[[2]]
  uhat <- as.complex(y[1:N] + (1i) * y[(N+1):length(y)])
  d_uhat <- -a^2 * (kappa^2) * uhat
  d_uhat_ri <- c(Re(d_uhat), Im(d_uhat))
  list(d_uhat_ri)
}

# basically function from np.fft.fftfreq
# Define discrete wavenumbers
kappa <- seq(-5, 5, length.out=(N+1))
kappa <- ifftshift(kappa[1:(length(kappa)-1)])
kappa <- kappa * 2 * pi


# Initial condition
u0 <- rep(0, length(x)-1)
u0[as.integer((L/2 - L/10)/dx):as.integer((L/2 + L/10)/dx)] <- 1

fft_ <- fft(u0) #Fourier transform of initial condition

# ode function doesn't play well with complex numbers, so we recast 
# the state u0hat from an N-element complex vector to a 2N-element real vector
u0hat_ri <- c(Re(fft_), Im(fft_))

# Simulate in Fourier frequency domain
dt <- 0.1
t <- seq(0,10,dt)

pars <- list()
pars[[1]] <- kappa
pars[[2]] <- a

uhat_ri <- ode(y = u0hat_ri, times = t, func = rhsHeat, parms = pars, method="ode45")

uhat <- matrix(rep(0, length(t)*N), nrow=length(t), ncol=N)

for (i in 2:length(t)) {
  uhat[i-1,] <- as.complex(uhat_ri[i,2:(N+1)] + (1i) * uhat_ri[(N+2):(2*N+1)])
}

u <- matrix(rep(0, length(t)*N), nrow=length(t), ncol=N)

for (i in 2:length(t)) {
  uhat[i-1,] <- as.complex(uhat_ri[i,2:(N+1)] + (1i) * uhat_ri[(N+2):(2*N+1)])
  u[i-1,] <- fft(uhat[i-1,], inverse=TRUE)/length(uhat[i-1,])
}

u <- Re(u)


u1 <- data.frame(x = 1, y = 1:length(u[1,]), cut = u[1,])

plot_ly(u1, x = ~x, y = ~y, z = ~cut, type = 'scatter3d', mode = 'lines', color = ~cut)

#dens <- with(diamonds, tapply(price, INDEX = cut, density))
#data <- data.frame(
#  x = unlist(lapply(dens, "[[", "x")),
#  y = unlist(lapply(dens, "[[", "y")),
#  cut = rep(names(dens), each = length(dens[[1]]$x)))

#plot_ly(data, x = ~x, y = ~y, z = ~cut, type = 'scatter3d', mode = 'lines', color = ~cut)
