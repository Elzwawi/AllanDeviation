# Session info
# sessionInfo()

# Libraries
library(ggplot2)
library(scales)

# reading and cleaning data
Phasenoise <- read.csv("PN.csv")        # reading file
Phasenoise <- Phasenoise[,1:2]                  # Cleaning data
colnames(Phasenoise) <- c("Freq","PN")

# Nominal Frequency
fo <- 10.949297E6

# Visualize data 
ggplot(Phasenoise, aes(x=Freq, y=PN)) + geom_line() +  scale_x_log10(breaks=10^(0:7)) + 
  xlab("Frequency (Hz)") + ylab("Phase noise (dBc/Hz)") + annotation_logticks(sides="b")

# Calculate Sy(f) from L(f)
Phasenoise$Lf <- 10^((Phasenoise$PN)/10)    # converting L(f) from dBz/Hz to normal scale
Phasenoise$SyF <- (2*(Phasenoise$Freq/fo)^2)*Phasenoise$Lf 
ggplot(Phasenoise, aes(x=Freq, y=SyF)) + geom_line() + scale_x_log10() + scale_y_log10() +
  xlab("Frequency (Hz)") + ylab("Freq. mod. noise (1/Hz)") + annotation_logticks(sides="bl")

# Calculate frequency intervals and mean Sy(F) in each interval
n <- nrow(Phasenoise)                # number of data points
df <- NULL                           # Variable for number of freq. intervals
SyFm <- NULL                         # Variable for mean Sy(F) ineach interval
m <- n-1

# df Array of frequency intervals (integration steps)
for (i in 1:n-1) {                  
 df[i] <- Phasenoise$Freq[i+1] - Phasenoise$Freq[i]
}


# mean fractional frequency fluctatuib Sy(F) in each integration interval
for (j in 1:n-1) {                  
  SyFm[j] <- (Phasenoise$SyF[j+1] + Phasenoise$SyF[j])/2
}

## Calculate integral (Allan deviation)
#tau <- c(0.0001, 0.001, 0.01, 0.1, 1, 5, 10, 50, 100, 1000)   # array of tau values
tau <- c(seq(from=0.0001, to=0.0009, by=0.0001), 
         seq(from=0.001, to=0.009, by=0.001),
         seq(from=0.01, to=0.09, by=0.01),
         seq(from=0.1, to=0.9, by=0.1),
         seq(from=1, to=9, by=1),
         seq(from=10, to=90, by=10),
         seq(from=100, to=900, by=100),
         seq(from=1000, to=9000, by=1000))

         # multiplier <- rbind(1:n, 1:n, 1:n, 1:n, 1:n)       # just defining multiplier to avoid error message (same length as tau)
multiplier <- matrix(1:(n*length(tau)),nrow=length(tau), ncol=n)

# calculate multiplier
# matrix with tau-rows, every row contains a vector corresponding to the measurement freq reange

for (j in 1:length(tau)) { 
  for (i in 1:n) {
    multiplier[j,i] <- ((sin(pi/180*pi*tau[j]*Phasenoise$Freq[i]))^4)/((pi*tau[j]*Phasenoise$Freq[i])^2)
  }
}

## Area under each interval
#integrand <- rbind(1:m, 1:m, 1:m, 1:m, 1:m)    # dummy place holder
integrand <- matrix(1:(m*length(tau)),nrow=length(tau), ncol=m)

for (j in 1:length(tau)) { 
  for (i in 1:m) {
integrand[j,i] <- multiplier[j,i]*SyFm[i]*df[i]
}
}

## Integration: Summing all areas
Allan_var <- 1:length(tau) # place holder for allan variance

# Calculating the integral for each value of tau
for (j in 1:length(tau)) { 
  Allan_var[j] <- 2*sum(integrand[j,])
}

# Allan deviation
Allan_dev <- sqrt(Allan_var)           

# Plot Allan deviation
data <- data.frame(cbind(tau,Allan_dev))
ggplot(data, aes(x=tau, y=Allan_dev)) + geom_line() + geom_point() + 
  scale_x_log10() + scale_y_log10() + xlab("tau (sec)") + ylab("Allan deviation") +
  annotation_logticks(sides="bl")

# writing a csv file
write.csv(data, file = "allan_dev.csv")




