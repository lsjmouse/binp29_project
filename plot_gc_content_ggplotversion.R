#prepare necessary package
library(seqinr)
library(ggplot2)

#load input file
ma<-read.fasta(file = "~/binp29project/microcystis aeruginosa.fasta")
maseq<-ma[[1]]
mv<-read.fasta(file="~/binp29project/Microcystis viridis.fasta")
mvseq<-mv[[1]]

#sliding window gc plot function
slidingwindowplot<-function(windowsize, inputseq, startpoints, endpoints){
  if (endpoints == 0) {
    starts<-seq(startpoints, length(inputseq)-windowsize, by=windowsize)
  }
  if (endpoints != 0) {
    starts<-seq(startpoints, endpoints-windowsize, by=windowsize)
  }
  n<-length(starts)
  chunkgcs<-rep(0,n)
  for (i in 1:n) {
    chunk<-inputseq[starts[i]:(starts[i]+windowsize-1)]
    chunkGC<-GC(chunk)
    print(chunkGC)
    chunkgcs[i]<-chunkGC
  }
  gc_data<-data.frame(starts,chunkgcs)
  plot<-ggplot(gc_data,aes(x=starts, y=chunkgcs)) +
    geom_point(colour="black",fill='green', pch=21, size=1) +
    geom_line(colour='red')+
    xlab("Nucleotied start position") +
    ylab("GC content")  +
    ggtitle("GC content plot") +
    theme_bw()
  print(plot)
}

#usage: slidingwindowplot (windowsize, inputseq, startpoints, endpoints)
#if you don't want to specific start position and end position, 
#just type 1 for start position and 0 for end position

#example
slidingwindowplot(50000,mvseq,1500000,0)
slidingwindowplot(50000,maseq,1500000,0)

