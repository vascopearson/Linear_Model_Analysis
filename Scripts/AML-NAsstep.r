library(ggplot2)
library(GGally)
library(car)
library(lmtest)
library(nortest)
library(calibrate)

#dados
life <- read.csv("https://web.tecnico.ulisboa.pt/~ist13493/AML2020_2021/Project/LifeExpectancyData.csv")
life <- as.data.frame(life)
life<-life[,-c(2,3)] # retirar variáveis "year" e "status"

#variável categórica
life$Country<-as.factor(life$Country)

#BMI < 40
life[,9]=ifelse(life[,9]>40,40,life[,9])

#Retirar NA's
life<-na.omit(life)
rownames(life)<-c(1:length(life[,1]))
summary(life)

#OUTLEIRS CARIA RAPIDO
out=boxplot(life)$out
out
life <- life[-which(life$Life.expectancy %in% out),]

#Reordenar
rownames(life)<-c(1:length(life[,1]))

#dados sem coluna dos países
life_sem_paises = life[,-1]

pairs(life)
ggpairs(life[,-1])+ scale_color_manual(values=c("#FA58AC"))

n<-dim(life)[1]

set.seed(3493)             # for reproducible example
test.ind<-sample(n,0.2*n)  # random sample of 20% of data

#Training data
data.train<-life[-test.ind,]
summary(data.train)
dim(data.train)

#Test data
data.test <- life[test.ind,]
summary(data.test)
dim(data.test)

#troca o nauru
#data.train<-rbind(data.train,data.test[478,])
#data.test<-data.test[-478,]

#Reordenar e retirar
rownames(data.test)<-c(1:length(data.test[,1]))
data.train<-rbind(data.train,data.test[54,]) #Ireland
data.test<-data.test[-54,]
rownames(data.test)<-c(1:length(data.test[,1]))
data.train<-rbind(data.train,data.test[54,]) #Equatorial Guine
data.test<-data.test[-54,]

#Anova
anova_alt = function (object, reg_collapse=TRUE,...) 
{
  if (length(list(object, ...)) > 1L) 
    return(anova.lmlist(object, ...))
  if (!inherits(object, "lm")) 
    warning("calling anova.lm(<fake-lm-object>) ...")
  w <- object$weights
  ssr <- sum(if (is.null(w)) object$residuals^2 else w * object$residuals^2)
  mss <- sum(if (is.null(w)) object$fitted.values^2 else w * 
               object$fitted.values^2)
  if (ssr < 1e-10 * mss) 
    warning("ANOVA F-tests on an essentially perfect fit are unreliable")
  dfr <- df.residual(object)
  p <- object$rank
  if (p > 0L) {
    p1 <- 1L:p
    comp <- object$effects[p1]
    asgn <- object$assign[stats:::qr.lm(object)$pivot][p1]
    nmeffects <- c("(Intercept)", attr(object$terms, "term.labels"))
    tlabels <- nmeffects[1 + unique(asgn)]
    ss <- c(vapply(split(comp^2, asgn), sum, 1), ssr)
    df <- c(lengths(split(asgn, asgn)), dfr)
    if(reg_collapse){
      if(attr(object$terms, "intercept")){
        collapse_p<-2:(length(ss)-1)
        ss<-c(ss[1],sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c(tlabels[1],"Source")
      } else{
        collapse_p<-1:(length(ss)-1)
        ss<-c(sum(ss[collapse_p]),ss[length(ss)])
        df<-c(df[1],sum(df[collapse_p]),df[length(df)])
        tlabels<-c("Regression")
      }
    }
  }else {
    ss <- ssr
    df <- dfr
    tlabels <- character()
    if(reg_collapse){
      collapse_p<-1:(length(ss)-1)
      ss<-c(sum(ss[collapse_p]),ss[length(ss)])
      df<-c(df[1],sum(df[collapse_p]),df[length(df)])
    }
  }
  
  ms <- ss/df
  f <- ms/(ssr/dfr)
  P <- pf(f, df, dfr, lower.tail = FALSE)
  table <- data.frame(df, ss, ms, f, P)
  table <- rbind(table, 
                 colSums(table))
  if (attr(object$terms, "intercept")){
    table$ss[nrow(table)]<- table$ss[nrow(table)] - table$ss[1]
  }
  table$ms[nrow(table)]<-table$ss[nrow(table)]/table$df[nrow(table)]
  table[length(P):(length(P)+1), 4:5] <- NA
  dimnames(table) <- list(c(tlabels, "Error","Total"), 
                          c("Df","SS", "MS", "F", 
                            "P"))
  if (attr(object$terms, "intercept")){
    table <- table[-1, ]
    table$MS[nrow(table)]<-table$MS[nrow(table)]*(table$Df[nrow(table)])/(table$Df[nrow(table)]-1)
    table$Df[nrow(table)]<-table$Df[nrow(table)]-1
  }
  structure(table, heading = c("Analysis of Variance Table\n"), 
            class = c("anova", "data.frame"))
}


m1 <- lm(data.train$Life.expectancy~.,data=data.train) # com todas as variaveis X
summary(m1)

m1_sempaises <- lm(data.train[,-1]$Life.expectancy~.,data=data.train[,-1]) # com todas as variaveis X


anova(m1_sempaises)
anova_alt(m1)

extractAIC(m1) #AIC
extractAIC(m1, k = log(n)) #BIC
sum((m1$residuals/( 1-hatvalues(m1)))^2) #PRESSp

#Previsão
y.pred.m1 <- predict(m1,data.test) 
summary(sqrt((y.pred.m1 - data.test$Life.expectancy)^2))

d<-data.frame(yp.m1=y.pred.m1, y=data.test$Life.expectancy)
ggplot(d, aes(yp.m1, y, color =yp.m1 )) +
  geom_point(shape = 16, size = 5, show.legend = FALSE) +
  theme_minimal()



######## step FOI FEITO SEM PAISES(ADDED AFTER)
m.base <-lm(data.train[,-1]$Life.expectancy~Adult.Mortality,data=data.train[,-1]) #É a variável inicial mais significativa
summary(m.base)

m.full<- lm(data.train[,-1]$Life.expectancy~.,data=data.train[,-1])  
summary(m.full)

step.forward.NA <- step(m.base,scope =list(upper=m.full,lower=~1), direction = "forward", trace=FALSE) #escolhe com base no AIC

step.backward.NA <- step(m.full, direction = "backward", trace=FALSE )

step.both.NA <- step(m.base, scope = list(upper=m.full, lower=~1 ), direction = "both", trace=FALSE)

#Backward elimination using p-values to delete predictors one-at-a-time
options(max.print=999999)
m.full<- lm(data.train[,-1]$Life.expectancy~.,data=data.train[,-1])
summary(m.full)

m.new = update(m.full, .~. - Measles  )
summary(m.new)

m.new = update(m.new, .~. - thinness..1.19.years  )
summary(m.new)

m.new = update(m.new, .~. - Population )
summary(m.new)

m.new = update(m.new, .~. - Alcohol  )
summary(m.new)

m.new = update(m.new, .~. - Hepatitis.B    )
summary(m.new)

m.new = update(m.new, .~. - Diphtheria    )
summary(m.new)

m.new = update(m.new, .~. - GDP    )
summary(m.new)

m.new = update(m.new, .~. - thinness.5.9.years    )
summary(m.new)

step.backward.p.NA<-m.new

library(MASS)
summary(m.full)
dropterm(m.full, test = "F" )

m.new<-update(m.full, .~. - Measles )
dropterm(m.new, test = "F" )

m.new<-update(m.new, .~. - thinness..1.19.years )
dropterm(m.new, test = "F" )

m.new<-update(m.new, .~. - Population )
dropterm(m.new, test = "F" )

m.new<-update(m.new, .~. - Alcohol  )
dropterm(m.new, test = "F" )

m.new<-update(m.new, .~. - Hepatitis.B  )
dropterm(m.new, test = "F" )

m.new<-update(m.new, .~. - Diphtheria   )
dropterm(m.new, test = "F" )

m.new<-update(m.new, .~. - GDP   )
dropterm(m.new, test = "F" )

m.new<-update(m.new, .~. - thinness.5.9.years   )
dropterm(m.new, test = "F" )

step.backward.f.NA <- m.new
## o step.backward.f.NA e o step.backward.p.NA sao o mesmo e a ordem é a mesma


library(MASS)
m.null<- lm(data.train[,-1]$Life.expectancy~ 1,data=data.train[,-1]) 
summary(m.null)

addterm(m.null, scope=m.full, test="F" )
m.new<-update(m.null, .~. + Schooling)
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + HIV.AIDS )
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + Adult.Mortality    )
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + Income.composition.of.resources)
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + percentage.expenditure)
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + BMI)
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + Polio)
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + Total.expenditure )
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + thinness..1.19.years )
addterm(m.new, scope=m.full, test="F" )

m.new<-update(m.new, .~. + Diphtheria   )
addterm(m.new, scope=m.full, test="F" )

step.forward.f.NA <- m.new

step.both.aic.NA <- step(m.base, direction = "both", scope= formula(m1_sempaises))

summary(step.backward.NA) #11 VAR BEST
anova_alt(step.backward.NA) # a
summary(step.backward.p.NA)#11 VAR pior #igual ao step backwards f
anova_alt(step.backward.p.NA) # b
summary(step.both.NA)#11 VAR pior
anova_alt(step.both.NA) # c
summary(step.both.aic.NA)#11 VAR pior
anova_alt(step.both.aic.NA) # d
summary(step.forward.NA)# 12 VAR pior
anova_alt(step.forward.NA) # e
summary(step.forward.f.NA)# 12 VAR pioe
anova_alt(step.forward.f.NA)# f 


best.fit<-c("Country", "Adult.Mortality", "percentage.expenditure", "BMI", "under.five.deaths", "Polio", "Total.expenditure", "HIV.AIDS", "thinness.5.9.years", "Income.composition.of.resources", "Schooling") 
#Retirámos o infant deaths pq ja tinhamos o  tinha mais erro
m.bestsubset <- lm(data.train$Life.expectancy ~ .,data=data.train[,best.fit])

summary(m.bestsubset)
anova_alt(m.bestsubset)

y.pred.m.bestsubset <- predict(m.bestsubset,data.test) 
summary(sqrt((y.pred.m.bestsubset - data.test$Life.expectancy)^2))

d<-data.frame(yp.m.bestsubset=y.pred.m.bestsubset, y=data.test$Life.expectancy)
ggplot(d, aes(yp.m.bestsubset, y, color =yp.m.bestsubset )) +
  geom_point(shape = 16, size = 5, show.legend = FALSE) +
  theme_minimal()

#Diagnostico
m.bestsubset <- lm(data.train$Life.expectancy ~ .,data=data.train[,best.fit]) # considere the best fit obtained by R^2 adjus criteria
plot(m.bestsubset) #diagnostic plots

#correlação entre erros
plot(residuals(m.bestsubset),ylab="Valor dos res?duos",xlab="?ndice",col=rgb(0.8,0.45,0.1),pch=16)
dwtest(m.bestsubset) #teste Durbin Watson
## FALHA

#Avaliar homoscedasticity 
ncvTest(m.bestsubset) # Breusch-Pagan test
### ACHO QUE A HOMOSCEDASTICITY NAO FALHA CARALHOOOO ### (Verifiquem)

#Avaliar normalidade dos erros e invariancia do valor esperado dos erros
boxplot(resid(m.bestsubset),col=rgb(0.5,0.8,1))
stem(fitted.values(m.bestsubset))
hist(resid(m.bestsubset),col=rgb(0.5,0.8,1),pch=16,main="",xlab="Res?duos",ylab="Frequ?ncia")
#Qui Quadrado
pearson.test(resid(m.bestsubset))

#shapiro teste para normalidade
shapiro.test(residuals(m.bestsubset))

vif(m.bestsubset)
crPlots(m.bestsubset)
library(gvlma)
gvmodel <- gvlma(m.bestsubset)
summary(gvmodel)


###########ELE NAO ESTA A DAR UPDATE E EU NAO SEI PORQUE
############## metodo para fazer com o do caria##############################
############## SO CORRER NA INICIALIZACAO:
newlife<-life
############## REPETIR ATE AGRADAR:
############## REPETIR APARTIR DA SEGUNDA ITERADA
newlife<- rbind(data.test,data.train)
#####################################################
n<-dim(newlife)[1]
rownames(newlife)<-c(1:length(newlife[,1]))
set.seed(3493)            
test.ind<-sample(n,0.2*n)  
data.train<-newlife[-test.ind,]
data.test <- newlife[test.ind,]
data.train<-rbind(data.train,data.test[478,])
data.test<-data.test[-478,]
rownames(data.train)<-c(1:length(data.train[,1]))
m.bestsubset <- lm(data.train$Life.expectancy ~ .,data=data.train[,best.fit]) 
extractAIC(m.bestsubset)#AIC
extractAIC(m.bestsubset, k = log(n))#BIC
sum((m.bestsubset$residuals/( 1-hatvalues(m.bestsubset)))^2)#PRESS 
y.pred.m.bestsubset <- predict(m.bestsubset,data.test) 
                      # rownames(data.test)<-c(1:length(data.test[,1]))
                      # data.train<-rbind(data.train,data.test[c(16),]) #Ireland
                      # data.test<-data.test[-c(16),]
# rownames(data.test)<-c(1:length(data.test[,1]))
# data.train<-rbind(data.train,data.test[207,]) #Equatorial Guine
# data.test<-data.test[-207,]
#Ver se tiramos com base no outliertest
outlierTest(m.bestsubset)
#plot(m.bestsubset)
coisa=outlierTest(m.bestsubset,n.max=3)
OUTLIERS=as.numeric(as.vector(names(coisa$bonf.p)))
data.train<-data.train[-c(OUTLIERS),]
OUTLIERS
#APOS CARIA RAPIDO 2� com unkowns e cook
# 38 183 141
# 128 326
# 128 326
# 131
# 848
#  643
# 173
# 126
#  615
#DDFITS
k <- length(m.bestsubset$coefficients)-1

#aproximação de 2*sqrt(k/n), com um increase no valor por não serem muitos dados
cv=2*sqrt(k/n)

plot(dffits(m.bestsubset), 
     ylab = "Standardized DFFITS", xlab = "Índice",col=rgb(0.3,0.7,1))
abline(h = cv, lty = 2)
abline(h = -cv, lty = 2)

textxy(as.numeric(names(dffits(m.bestsubset)[which(dffits(m.bestsubset) < -cv | dffits(m.bestsubset) > cv)])), 
       dffits(m.bestsubset)[which(dffits(m.bestsubset) < -cv | dffits(m.bestsubset) > cv)], 
       as.numeric(names(dffits(m.bestsubset)[which(dffits(m.bestsubset) < -cv | dffits(m.bestsubset) > cv)])), cex=0.8,offset = -1)

#Retirar DFFITS 1�
DFFITS=which(abs(dffits(m.bestsubset))>cv)
DFFITS
# 42   77   78   80   85   87  101  108  116  121  122  169  171  173  178  179  189  190  234 
# 35   62   63   64   69   70   81   87   95  100  101  136  137  138  141  142  150  151  189 
# 240  255  353  354  369  376  379  385  391  392  393  399  421  444  518  519  525  542  545 
# 194  206  287  288  299  304  307  311  316  317  318  322  343  364  428  429  433  446  448 
# 592  599  607  629  631  636  657  660  670  673  681  690  692  694  713  715  716  723  725 
# 485  488  492  513  515  519  533  536  543  546  554  560  561  562  574  575  576  582  583 
# 736  819  837  841  850  851  854  858  861  920  948  969  975  980  981  982  986  987 1008 
# 591  654  668  671  677  678  681  682  685  735  757  776  780  782  783  784  787  788  805 
# 1016 1037 1054 1055 1078 1092 1093 1095 1097 1098 1110 1111 1209 1210 1211 1216 1217 1220 1263 
# 812  830  845  846  861  870  871  872  873  874  884  885  968  969  970  975  976  978 1010 
# 1268 1269 1270 1274 
# 1014 1015 1016 1020 

data.train<-data.train[-c(DFFITS),]

leveragePlots(m.bestsubset)# leverage plots 

# Cook's D plot
# identify D values > 4/(n-k-1)
cutoff <- 4/((nrow(data.train)-length(m.bestsubset$coefficients)-2))
plot(m.bestsubset, which=4, cook.levels=cutoff)
abline(h = cutoff, lty = 2)

# Influence Plot
influencePlot(m.bestsubset, main="Influence Plot", sub="Circle size is proportial to Cook's Distance")
#tirando cook e unkowns
data.train<-data.train[-c(141,183,180,49,69),]
data.train<-data.train[-c(102,793,326,142,224),]
data.train<-data.train[-c(197,267,378,185),]
data.train<-data.train[-c(259,628,270,931,932),]
data.train<-data.train[-c(299,175,283,860,917),]
data.train<-data.train[-c(364,471,501,717,924),]
data.train<-data.train[-c(131,219,384,644,864),]
data.train<-data.train[-c(77,120,728),]
data.train<-data.train[-c(271,376,624,912,310),]
data.train<-data.train[-c(154,287,415,916),]
data.train<-data.train[-c(254,550,769,907,910),]
data.train<-data.train[-c(178,602,795,53,174),]
data.train<-data.train[-c(142,725,789,84,132),]
data.train<-data.train[-c(237,471,556,842,901),]
data.train<-data.train[-c(328,616,642,106,111),]
data.train<-data.train[-c(14,173,622,447,659),]
data.train<-data.train[-c(660,790,824,889,700),]
data.train<-data.train[-c(135,488,211,104,891),]
data.train<-data.train[-c(831,553,104,23,831),]
data.train<-data.train[-c(203,126,239),]
data.train<-data.train[-c(522,595,814,240),]
data.train<-data.train[-c(100,466,728,163,440),]
data.train<-data.train[-c(68,138,465),]
data.train<-data.train[-c(152,460,590),]
data.train<-data.train[-c(455,824,473,11,341),]
data.train<-data.train[-c(281,71,59,32,267),]
data.train<-data.train[-c(237,292,720),]
data.train<-data.train[-c(323,448,835,196),]
data.train<-data.train[-c(596,9,564),]
data.train<-data.train[-c(15,716,820),]
data.train<-data.train[-c(80,115,172,83,409),]
data.train<-data.train[-c(2,752,755,165),]
data.train<-data.train[-c(113,173,272,592,676),]
data.train<-data.train[-c(126,93,166,2,205),]
data.train<-data.train[-c(694,177,128,420,293),]
data.train<-data.train[-c(210,249,326),]
data.train<-data.train[-c(26,135),]
data.train<-data.train[-c(126,586,396),]
data.train<-data.train[-c(224,453,701,66,831),]
data.train<-data.train[-c(74,83,350,54,675),]
data.train<-data.train[-c(155,584,812,821),]
data.train<-data.train[-c(131,513,689),]
data.train<-data.train[-c(186,231,414,85,568),]
data.train<-data.train[-c(160,106,750,158,161),]
data.train<-data.train[-c(22,306,469),]
data.train<-data.train[-c(155,154,691,813),]
data.train<-data.train[-c(389,490,623,811),]
data.train<-data.train[-c(127,157,740,111),]
data.train<-data.train[-c(388,167,291,96),]
data.train<-data.train[-c(217,460,615,728),]
data.train<-data.train[-c(32,90,726,181),]
data.train<-data.train[-c(72,521,575),]
data.train<-data.train[-c(521),]
data.train<-data.train[-c(314,350,523),]
data.train<-data.train[-c(120,136),]
data.train<-data.train[-c(82,488,561),]
data.train<-data.train[-c(80,169,255,608,790),]
data.train<-data.train[-c(96,24,114,200),]
data.train<-data.train[-c(415,538,712,189,477),]
data.train<-data.train[-c(36,433,558),]
data.train<-data.train[-c(159),]
data.train<-data.train[-c(27,284,391,248,780),]
data.train<-data.train[-c(56,78,527),]
data.train<-data.train[-c(314,221,392),]
data.train<-data.train[-c(249,501,402,773),]
data.train<-data.train[-c(12,657,152),]
data.train<-data.train[-c(540,476,163,168),]
data.train<-data.train[-c(106,649,749),]
data.train<-data.train[-c(7,304,391,765),]
data.train<-data.train[-c(626,471,489,676,765),]
data.train<-data.train[-c(626,471,489,676,765),]
data.train<-data.train[-c(41,311,421,653),]
data.train<-data.train[-c(66,536,677,764),]
data.train<-data.train[-c(706,582,635,508),]
data.train<-data.train[-c(135,192,450),]
data.train<-data.train[-c(258,340,358),]
data.train<-data.train[-c(207,509,682),]
data.train<-data.train[-c(446,535,556),]
data.train<-data.train[-c(229,381,444),]
data.train<-data.train[-c(270,567,651,602),]
data.train<-data.train[-c(57,141,297),]
data.train<-data.train[-c(405,707,249,94),]
data.train<-data.train[-c(748,134,201,404),]
data.train<-data.train[-c(37,127,676,660),]
data.train<-data.train[-c(117,504,611,671),]
data.train<-data.train[-c(154,264,363),]
data.train<-data.train[-c(685,621,136),]
data.train<-data.train[-c(544,362,487,553),]
data.train<-data.train[-c(50,100,231),]
data.train<-data.train[-c(512,232,534,733),]
data.train<-data.train[-c(255,209,291),]
data.train<-data.train[-c(182,553,665,271),]
data.train<-data.train[-c(30,177,578,724),]
data.train<-data.train[-c(26,586,674),]
data.train<-data.train[-c(257,617,625),]
data.train<-data.train[-c(499,2,71,286),]
data.train<-data.train[-c(239,294,414),]
data.train<-data.train[-c(326,445,681),]
data.train<-data.train[-c(686,719,65,515),]
data.train<-data.train[-c(111,376,373),]
data.train<-data.train[-c(19,114,180),]
data.train<-data.train[-c(266,445,524,683),]
data.train<-data.train[-c(214,377,441,151),]
data.train<-data.train[-c(65,461,540),]
data.train<-data.train[-c(558,20,213,277),]
data.train<-data.train[-c(69,89,636),]
data.train<-data.train[-c(41,2,241),]
data.train<-data.train[-c(276,60,79,168),]
data.train<-data.train[-c(28,329,131),]
data.train<-data.train[-c(699,182,301,460),]
data.train<-data.train[-c(123,225,280,417),]
data.train<-data.train[-c(43,131,620),]
data.train<-data.train[-c(156,430),]
data.train<-data.train[-c(113,24,483,606),]
data.train<-data.train[-c(234,125,596),]
data.train<-data.train[-c(104,194,215),]
data.train<-data.train[-c(685,113,28,675),]
data.train<-data.train[-c(573,119,242,495),]
data.train<-data.train[-c(27,11,130,280),]

lifesemout<-newlife