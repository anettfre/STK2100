data = read.csv('leukemia_big', header = TRUE)
datat = t(data) #transponerer siden radene er observasjonene

set.seed(2100)

testnames = c("ALL.4", "ALL.8", "ALL.10", "ALL.11", "ALL.13", "ALL.18", "AML", "AML.1",
  "AML.4", "AML.6", "AML.8", "ALL.23", "ALL.26", "ALL.29", "ALL.31", "ALL.32", "ALL.35",
  "ALL.39", "ALL.40", "ALL.41", "ALL.42", "AML.16", "AML.22", "AML.24")

test_true = is.element(rownames(datat), testnames)
test.data = datat[test_true,]
train.data = datat[!test_true,]

data.type = grepl("ALL", rownames(datat))
test.y = data.type[test_true]
train.y = data.type[!test_true]
#a.1
library(glmnet)
lambdas = exp(-7:2)
lasso =  cv.glmnet(x = train.data, y = train.y,  nfolds = 48, alpha = 1, standardize = TRUE, family = "binomial", lambda = lambdas, grouped = FALSE)
lasso$lambda.min

mod.ridge = cv.glmnet(x = train.data, y = train.y, nfolds = 48, alpha = 0, standardize = TRUE, family = "binomial", lambda = lambdas, grouped = FALSE)
mod.ridge$lambda.min

#a.2
#lasso
'''
lasso.predict.train = predict(lasso$glmnet.fit, newx = train.data, s = lasso$lambda.min, type= "response")
lasso.hat.train = lasso.predict.train >= 0.5

lasso.predict.test = predict(lasso$glmnet.fit, newx = test.data, s = lasso$lambda.min, type= "response")
lasso.hat.test = lasso.predict.test >= 0.5
                        
error.lasso.train = mean((train.y - lasso.hat.train)^2)
error.lasso.test = mean((test.y - lasso.hat.test)^2)


#ridge
ridge.predict.train = predict(mod.ridge$glmnet.fit, newx = train.data, s = mod.ridge$lambda.min,  type= "response")
ridge.hat.train = ridge.predict.train >= 0.5

ridge.predict.test = predict(mod.ridge$glmnet.fit, newx = test.data, s = mod.ridge$lambda.min,  type="response")
ridge.hat.test = ridge.predict.test >= 0.5

error.ridge.train = mean((train.y - ridge.hat.train)^2)
error.ridge.test = mean((test.y - ridge.hat.test)^2)

cbind(error.lasso.train, error.lasso.test, error.ridge.train, error.lasso.test)
'''

mse.logistic.cv = function(cv.mod, x, y) {
  y.probs = predict(cv.mod$glmnet.fit, s=cv.mod$lambda.min, newx=x, type="response")
  y.hat = y.probs >= 0.5
  return(mean((y - y.hat)^2))
}

cbind(mse.logistic.cv(lasso, train.data,train.y), mse.logistic.cv(lasso, test.data, test.y), mse.logistic.cv(mod.ridge, train.data,train.y), mse.logistic.cv(mod.ridge, test.data, test.y))


#c
ALL = grepl("ALL", rownames((datat)))
AML = grepl("AML", rownames((datat)))
p.values <- vector(mode='numeric', length=ncol(datat))
for (i in (1:ncol(datat))) {
  t.test.result = t.test(datat[ALL, i], datat[AML, i])
  p.values[i] = t.test.result$p.value
}
best.p.values = order(p.values)[1:9]
p.values[best.p.values]

test.c = test.data[, best.p.values]
train.c = train.data[, best.p.values]
log.model = glm(train.y ~ ., data = as.data.frame(train.c), family = 'binomial')

log.predict.train = predict(log.model, newdata = as.data.frame(train.c),  type="response")
log.hat.train = log.predict.train >= 0.5
error.logistic.train = mean((train.y - log.hat.train)^2)

log.predict.test = predict(log.model, newdata = as.data.frame(test.c),  type="response")
log.hat.test = log.predict.test >= 0.5
error.logistic.test = mean((test.y - log.hat.test)^2)

cbind(error.logistic.train, error.logistic.test)


#c.3
ALL2 = grepl("ALL", rownames((train.data)))
AML2 = grepl("AML", rownames((train.data)))
p.values2 <- vector(mode='numeric', length=ncol(train.data))
for (i in (1:ncol(train.data))) {
  t.test.result = t.test(train.data[ALL2, i], train.data[AML2, i])
  p.values2[i] = t.test.result$p.value
}
best.p.values2 = order(p.values2)[1:9]
best.p.values2
p.values2[best.p.values2]

test.c2 = test.data[, best.p.values2]
train.c2 = train.data[, best.p.values2]
log.model2 = glm(train.y ~ ., data = as.data.frame(train.c2), family = 'binomial')

log.predict.train2 = predict(log.model2, newdata = as.data.frame(train.c2),  type="response")
log.hat.train2 = log.predict.train2 >= 0.5
error.logistic.train2 = mean((train.y - log.hat.train2)^2)

log.predict.test2 = predict(log.model2, newdata = as.data.frame(test.c2),  type="response")
log.hat.test2 = log.predict.test2 >= 0.5
error.logistic.test2 = mean((test.y - log.hat.test2)^2)

cbind(error.logistic.train2, error.logistic.test2)


#d
init.observation = c("ALL", "ALL.10", "ALL.11", "AML", "AML.2")
init = which(is.element(rownames(datat), init.observation))
for (i in 2:length(init)) {
  for (j in 1:(i-1)) {
    place = init[c(i, j)]
    k1 = kmeans(datat, datat[place,])
    mod.k1 = glm(ALL ~ k1$cluster, family = 'binomial')
    print(i)
    print(j)
    print(mean((ALL - (mod.k1$fitted.values >= 0.5))^2))
  }
}

#e
d = dist(datat)
h1 = hclust(d, method = "complete")
plot(h1, xlab = "", ylab = "", axes = FALSE, main = "", sub = "", lwd = 2)
#plot(h1, cex = 0.6, hang = -1)
#rect.hclust(h1, k = 5, border = 2:6)


pc = prcomp(datat, center = TRUE, scale = TRUE)
d2 = dist(pc$x[,1:2])
h2 = hclust(d2, method = "complete")
plot(h2, xlab = "", ylab = "", axes = FALSE, main = "", sub = "", lwd = 2)
plot(h2, cex = 0.6, hang = -1)
#rect.hclust(h2, k = 2, border = 2:4)

library(factoextra)

fviz_eig(pc, ncp = 72)
colDen <- ifelse(data.type, 1, 2)
ColorDendrogram(h1, colDen, main = "", branchlength = 7.5, sub="cluster e.1",ylab = "")
legend("topright", legend=c("ALL", "AML"),
       col=c("red", "green"), lty=1, cex=0.5)
ColorDendrogram(h2, colDen, main = "", branchlength = 21, sub="cluste e.3r",ylab = "")
legend("topright", legend=c("ALL", "AML"),
       col=c("red", "green"), lty=1, cex=0.5)
