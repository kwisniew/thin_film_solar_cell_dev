read.table("IV_data.txt", header = FALSE, sep = "", dec = ".") -> IV
IV_curve = data.frame(voltage = IV$V1, current = IV$V2)
IV_curve = IV_curve[order(IV_curve$voltage),]
IV_curve$voltage = IV_curve$voltage-IV_curve$voltage[which.min(IV_curve$current)]
plot(IV_curve, log="y")

read.table("CV_data.txt", header = FALSE, sep = "", dec = ".") -> CV
CV_curve = data.frame(voltage = CV$V1, capacitance = CV$V2)
CV_curve = CV_curve[order(CV_curve$voltage),]
CV_curve$capacitance = abs(CV_curve$capacitance)
CV_curve$voltage = CV_curve$voltage-IV_curve$voltage[which.min(IV_curve$current)]
plot(CV_curve)


IV_curve$current[1:(which.min(IV_curve$current)-1)] = -IV_curve$current[1:(which.min(IV_curve$current)-1)]
plot(IV_curve)

plot(-CV_curve$voltage[1:11],(1/CV_curve$capacitance^2)[1:11],pch=19)
x <- -CV_curve$voltage[1:11]
y <- (1/CV_curve$capacitance^2)[1:11]
fit  <- lm(y~x)
seq(-0.7,0.7,0.1) -> xx
lines(xx, predict(fit, data.frame(x=xx)), col="red")
build_in_potential <- fit$coefficients[1]/fit$coefficients[2]

eps0 <- 8.85*10^-14
epsRelativ <- 13.0
S<-0.5*0.5
q<-1.602176565*10^-19
Na<- 10^(14)
Nd<- 10^(14)
N=(1/Na+1/Nd)^-1
N_derived = 2/(fit$coefficients[2]*q*epsRelativ*eps0*S^2)