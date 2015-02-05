args<-commandArgs(trailingOnly=TRUE)
data<-read.table("undecided_reads", header=FALSE, sep="\t")

len<-length(data[,2])
for (x in 1:len) {
    lh_sv<-dnorm(data[x,2], mean=as.numeric(args[1]), sd=as.numeric(args[2]), log=TRUE) + ppois(data[x,3], 1, lower.tail=FALSE, log.p=TRUE)
    lh_ref<-dnorm(data[x,4], mean=as.numeric(args[1]), sd=as.numeric(args[2]), log=TRUE) + ppois(data[x,5], 1, lower.tail=FALSE, log.p=TRUE)
    if (abs(lh_sv-lh_ref) <= 0.6) {write(as.character(data[x,1]), file="still_undecided.txt", append=TRUE)}
    else if (lh_sv < lh_ref) {write(as.character(data[x,1]), file="reads_change.txt", append=TRUE)}
}
