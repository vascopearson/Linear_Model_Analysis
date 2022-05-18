library(IMTest)
library(nortest)
library(QuantPsyc)
library(jtools)
library(olsrr)

##########################  TESTES PARA AS SUPOSIÇÕES ##########################


pearson.test(resid(m.bestsubset)) #Qui Quadrado


shapiro.test(residuals(m.bestsubset)) #shapiro teste para normalidade


durbinWatsonTest(m.bestsubset) #teste Durbin Watson


ncvTest(m.bestsubset) #Avaliar homoscedasticity com Breusch-Pagan test


##############################  ANALISAR O MODELO ##############################
anova_alt(m.bestsubset) #MSE e +

cena=summary(m.bestsubset)
c(cena$adj.r.squared,cena$r.squared) #R^2 ajustado e R^2

extractAIC(m.bestsubset)#AIC

extractAIC(m.bestsubset, k = log(n))#BIC

sum((m.bestsubset$residuals/( 1-hatvalues(m.bestsubset)))^2)#PRESS (pode dar inf)?

ols_mallows_cp(m.bestsubset,m.full)#Cp

# previsão
previsao.y <- predict(m.bestsubset,data.test)
summary(sqrt((previsao.y - data.test$Life.expectancy)^2))

d<-data.frame(previsao.y, y=data.test$Life.expectancy) 
ggplot(d, aes(previsao.y, y, color =previsao.y )) +
  geom_point(shape = 16, size = 5, show.legend = FALSE) +
  theme_minimal()




vif(m.bestsubset)





m.prediction<-lm(data.test$Life.expectancy~.,data=data.test[best.fitcp])




previsao.y <- predict(m.bestsubset,data.test)
summary(sqrt((previsao.y - data.test$Life.expectancy)^2)) #se min e max forem mto afastado o modelo n presta pq o valor previsto e observado para uma determinada linha vao ser mto diferentes. Nota que o R^2 era ''50%'', ou seja o m.full só explica 50% da variabilidade do y(life expectancy).

d<-data.frame(previsao.y, y=data.test$Life.expectancy) #gráfico que faz y(observados) em função de yp(previstos). Se os dados tiverem sob uma reta y=x tava bom mas nao vai tar
ggplot(d, aes(previsao.y, y, color =previsao.y )) +
  geom_point(shape = 16, size = 5, show.legend = FALSE) +
  theme_minimal()














