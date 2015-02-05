d<-scan(file="distance", what=double(), sep="\n")

qs<-quantile(d, c(0.01,0.99))
d1<-d[d > qs[1] & d < qs[2]]

op<-mat.or.vec(1,2)
op[1]<-mean(d1)
op[2]<-sd(d1)
op[1]
op[2]

write(op, file="dist_mean_sd.txt", sep="\n")
