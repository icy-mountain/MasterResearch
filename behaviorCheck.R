# model(freq=freq2), detail(model="LSI1")
# i=1 j=1
exp(0.03044 + 1.03211 + 1.03211 + (1-1)*(-0.39))
# i=4 j=4
exp(0.03044 + (4-4)*(-0.39))
model_foo <- globalAnalysResults[["LSI1"]]
summary(model_foo)
aic_foo <- model_foo$aic
lnL_foo <- logLik(model_foo)
print(sprintf("logLik:%s -2*logLik:%s AIC-(-2*logLik)=2k?:%s",lnL_foo, -2*lnL_foo, aic_foo+2*lnL_foo))
print(length(model_foo$coefficients))# k しかし係数がNAになってないかチェックして
print(model_foo$coefficients)
# quasi check
freq_q <- c(6, 2, 3, 1,
           9, 9, 2, 1,
           9, 2, 6, 1,
           12, 1, 2, 9)

print(model_foo$deviance)
print(model_foo$df.residual)
print(length(model_foo$coefficients))
