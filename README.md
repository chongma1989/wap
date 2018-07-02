# R package wap: Weight Adjusted P-value 

The R package wap implements a weight adjusted p-value approach to increase the power while controlling the familywise error rate or false discovery rate. 

## Install wap
* Download the zip file `wap_1.0.0.tar.gz`.
* In R console, run `install.packages("wap_1.0.0.tar.gz",repos=NULL,type="source")`. 

## Example 
```
X=matrix(rnorm(500),50)
beta1=0.5;beta2=0.1
Z1=X[,1]*beta1+rnorm(50)
Y=Z1*beta2+rnorm(50)
Z=cbind(Z1,matrix(rnorm(50000),50))
wpvalue(Y,X,Z,method="max",ncores=4)

```

