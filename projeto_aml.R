install.packages("gvlma")
install.packages("MASS")
install.packages("car")
install.packages("ggplot2")
install.packages("carData")
install.packages("VGAM")


library(MASS)
library(gvlma)
library(ggplot2)
library(car)
library(VGAM)

#2 TRATAMENTO DE DADOS


#2.1 Tratamento Preliminar dos Dados

d <- read.table("Trabalho2_Dados.txt") 
d[,8] <- as.factor(d[,8])
d[,9] <- as.factor(d[,9])
y <- d$V4
dados <- data.frame(y, x1=d$V1, x2=d$V2, x3=d$V3, x4=d$V5, x5=d$V6, x6=d$V7, x7=d$V8, x8=d$V9, x9=d$V10, x10=d$V11, x11=d$V12)


# Conjunto de Treino e Conjunto de Teste
s <- sample(113,22)
S = c(6, 87, 3, 35, 25, 97, 49, 59, 108, 70, 99, 76, 17, 61, 106,79, 9, 96, 14, 85, 113, 8)
teste.dados <- dados[S,]
treino.dados <- dados[-S,]
y <- treino.dados$y


# Modelo Completo
n<-dim(treino.dados)[1]
modc=lm(y ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11, treino.dados)
gvlma(modc)
 # suposições 
plot(modc) # 



#2.2 Verificar Normalidade da Variável Resposta

# Histograma
ggplot(treino.dados, aes(x=y)) + geom_histogram(color = "black", fill = "slategray4") + theme_minimal()

# Teste Shapiro-Wilk - antes de transformação
shapiro.test(treino.dados$y)


# sin(y) 
y1 <- sin(y)
qqnorm(y1, main = "Normal Q-Q",
       xlab = "Standardized Residuals", ylab = "Theoretical Quantiles")
qqline(y1, datax = FALSE, distribution = qnorm,probs = c(0.25, 0.75), lty = 2, qtype = 7)
ggplot(treino.dados, aes(x=y1)) + geom_histogram(color = "black", fill = "slategray4") + theme_minimal()

# Transformação Box-Cox
pt <-powerTransform(treino.dados$y, family = "bcPower")
lambda <- pt$lambda
yf <- treino.dados$y^lambda
#ggplot(treino.dados, aes(x=yf)) + geom_histogram(color = "black", fill = "slategray4")+ theme_bw()

# Teste Shapiro-Wilk - depois de transformação
shapiro.test(yf)


# Teste F de Ajustamento (Lack of Fit) - em substutuição do de cima




#3 MELHOR MODELO DE REGRESSÃO LINEAR

#3.1 Métodos de Eliminação de Covariáveis

# Novos Dados
dadosf <- data.frame(yf, x1=treino.dados$x1, x2=treino.dados$x2, x3=treino.dados$x3, x4=treino.dados$x4, x5=treino.dados$x5, x6=treino.dados$x6, x7=treino.dados$x7, x8=treino.dados$x8, x9=treino.dados$x9, x10=treino.dados$x10, x11=treino.dados$x11)
modcf=lm(yf ~ x1 + x2 + x3 + x4 + x5 + x6 + x7 + x8 + x9 + x10 + x11, dadosf)
summary(modcf)
anova(modcf)
extractAIC(modcf) # criterio AICp
extractAIC(modcf, k = log(n)) # criterio BICp
plot(modcf)
gvlma(modcf) 


# Métodos

#forward
fit <- lm(yf ~ 1, dadosf) #modelo base
modfow <- stepAIC(fit, direction ="forward", scope=list(upper=modcf, lower=fit)) 
summary(modfow)

#backward
modback <- stepAIC(modcf, direction = "backward") 
summary(modback)

#stepwise
modstep <- step(modcf)
summary(modstep)

# Escolhido
modr <- modfow

#Teste de Hipóteses: F parcial
anova(modr,modcf) 

#info modelo completo
Xc <- model.matrix(modcf)
Hc <- Xc%*%solve(t(Xc)%*%Xc)%*%t(Xc)
ssec <- sum((fitted(modcf)- yf)^2)
ssrc <- sum((fitted(modcf) - mean(yf))^2)

#info modelo reduzido
anova(modr)
summary(modr)
extractAIC(modr)
extractAIC(modr, k = log(n))
Xr <- model.matrix(modr)
Hr <- Xr%*%solve(t(Xr)%*%Xr)%*%t(Xr)
sser <- sum((fitted(modr)- yf)^2)
ssrr <- sum((fitted(modr) - mean(yf))^2)






#3.2 Variáveis de Segunda Ordem

dados2 <- data.frame(yf, x2=treino.dados$x2, x3=treino.dados$x3, x4=treino.dados$x4, x5=treino.dados$x5, x8=treino.dados$x8, x11=treino.dados$x11)
modc2 <- lm(yf ~(.)^2, dados2)

# Métodos

modfow2 <- stepAIC(modr, direction ="forward", scope=list(upper=modc2, lower=modr)) ##ESCOLHIDO
summary(modfow2)
extractAIC(modfow2) 
extractAIC(modfow2, k = log(n)) 

modback2 <- stepAIC(modc2, direction ="backward")
summary(modback2)
extractAIC(modback2) 
extractAIC(modback2, k = log(n)) 

modstep2 <- step(modc2)
summary(modstep2)
extractAIC(modstep2) 
extractAIC(modstep2, k = log(n)) 

# Modelo Final (com outliers)
modr2 <- modfow2



#3.4 Outliers e Pontos Influentes

# info modelo final
ssr <- sum((fitted(modr2) - mean(yf))^2)
sse <- sum((fitted(modr2) - yf)^2)
sst = ssr +sse 
mse= sse/tabela$Df[length(tabela$Df)]
msr = ssr/sum(tabela$Df[-length(tabela$Df)])

# Visualizar possíveis Outliers
boxplot(yf, main="Boxplot")

# Pontos Influentes
influencePlot(modr2, fill.col=carPalette()[1]) ##########
which(apply(influence.measures(modr2)$is.inf, 1, any))


# Alguns Resíduos

# habituais
X <- model.matrix(modr2)
H <- X%*%solve(t(X)%*%X)%*%t(X)
e <- yf - H%*%yf

# matriz de covariância dos resíduos
rcovmatrix <- mse * (diag(n) - H)
s2 <- diag(rcovmatrix)

# standerdized residual
stanr <- e/sqrt(mse)
which(sapply(stanr, function(x) x^2 > 3*sqrt(mse)))

# internally studentized residual
SR <- sqrt(s2)
istudr <- e/SR
denom <- diag(diag(n)-H)
d <- e/denom 
which(sapply(istudr-d, function(x) abs(x) > 1)) 

# hampel filter
mediana <- median(yf)
MAD <- mad(yf, center=median(yf))
upper.bound <- mediana + 3*MAD
lower.bound <- mediana - 3*MAD
upper.bound
lower.bound
which(sapply(yf, function(x) x>upper.bound || x< lower.bound)) 


# Retirar Outliers 

dadosout <- dadosf[-c(14,44,45,84),] #PORQUÊ O 29???

# modelo obtido anteriormente
modout <- lm(yf ~ x4 + x2 + x11 + x8 + x5 + x3 + x4:x11 + x11:x8 + x5:x3, data = dadosout)
anova(modout)
ssrout <- sum((fitted(modout) - mean(yf))^2)
summary(modout)
plot(modout)
extractAIC(modout) 
extractAIC(modout, k = log(n)) 


# modelo novo com dados sem outliers
novo.modout <- lm(yf ~. , dadosout)

#Eliminação de Covariáveis

#primeira ordem
fitout <- lm(yf ~ 1, dadosout) 
modfowout <- stepAIC(fitout, direction ="forward", scope=list(upper=novo.modout, lower=fitout)) 
summary(modfowout)
modrout <- modfowout

#segunda ordem
dadosout2 <- data.frame(yout=dadosout$yf, x2=dadosout$x2, x3=dadosout$x3, x4=dadosout$x4, x5=dadosout$x5, x8=dadosout$x8, x11=dadosout$x11)
modout2 <- lm(yout~(.)^2, dadosout2)
modfowout2 <- stepAIC(modrout, direction ="forward", scope=list(upper=modout2, lower=modrout)) ##ESCOLHIDO
summary(modfowout2)
extractAIC(modfowout2) 
extractAIC(modfowout2, k = log(n)) 
final.mod <- modfowout2




#4 TÉCNICAS DE DIAGNÓSTICO

#4.1 Validação do Modelo

#################################
################################
################################
#################################
################################
sink(file = "gvlmafinal_output.txt")
gvlma(final.mod) 
sink(file = NULL)
anova(final.mod)
summary(final.mod)

# Variance Inflation Factor - teste de multicolinearidade
v <-vif(final.mod, type="predictor") 
mean(v$`GVIF^(1/(2*Df))`) 

# Homocedasticidade dos Resíduos - teste Breusch-Pagan
ncvTest(final.mod) 

# Autocorrelação dos Resíduos - teste Durbin-Watson
durbinWatsonTest(final.mod) 


#4.2 Capacidade de Previsão

# RMSPE
expected <-predict(final.mod, newdata= teste.dados)
expected[7] <- 0
expected2 = expected^(1/lambda)
yteste <- teste.dados$y
SSPE <- (expected2- yteste)^2
MSPE <- sum(SSPE)/length(yteste)
RMSPE <-sqrt(MSPE)

# Range de Observações Originais
range(dados$y)

# Gráfico de Previsão
teste <- data.frame(real=yteste, expected=expected2)
ggplot(teste, aes(expected, real, color=real),scale) +
  geom_point(color="slategray4" ,shape = 16, size = 4,show.legend = FALSE) +
  coord_cartesian(x=c(0,8),y=c(0,8))+
  theme_minimal() + labs(title = "Prediction Best Model")+
  geom_abline(slope=1, intercept=0, color="grey10")



