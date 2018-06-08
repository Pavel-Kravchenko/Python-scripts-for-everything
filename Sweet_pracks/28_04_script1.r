attach(airquality)
x1 = Temp[Month==7]
x2 = Temp[Month==8]
r = cor(x1,x2)
cor.test(x1,x2, alternative = 'two.sided', method='pearson', conf.level = 0.95)
plot(x1,x2)
fischer <- function(r){0.5*log((1+r)/(1-r))}
z = qnorm(1-0.05/2)
n = 31
z.u = fischer(r) + z/sqrt(n-3)
z.l = fischer(r) - z/sqrt(n-3)

ht <- function(z){(exp(2*z)-1)/(exp(2*z)+1)}
ht(z.u)
ht(z.l)

cor.test(x1,x2, alternative = 'two.sided', method='spearman', conf.level = 0.95)

x = rnorm(100)
y = log(abs(x)) + rnorm(100, mean=0, sd = 0.5)
x1 = x[x>0]
y1 = y[x>0]
plot(x1,y1)
cor.test(x1,y1, alternative = 'two.sided', method='pearson', conf.level = 0.95)
cor.test(x1,y1, alternative = 'two.sided', method='spearman', conf.level = 0.95)
cor.test(x1,y1, alternative = 'two.sided', method='kendall', conf.level = 0.95)

detach(airquality)
