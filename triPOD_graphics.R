#!/usr/bin/Rscript

############################################################################### 
#
# The triPOD software detects chromosomal abnormalities using 
# the Parent-of-Origin-based Detection (POD) method.
# Author and maintainer: Joseph D. Baugher. 
# Email:jbaughe2(at)jhmi.edu.
#
# This file is part of the triPOD software.
# triPOD is free software: you can redistribute it and/or modify
#	it under the terms of the GNU General Public License as published by
#	the Free Software Foundation, either version 3 of the License, or
#	(at your option) any later version.
# triPOD is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
# You should have received a copy of the GNU General Public License
#   along with triPOD.  If not, see <http://www.gnu.org/licenses/>.
#
###############################################################################


library(shape)
library(TTR)

# Get arguments passed from Perl script
args           <- commandArgs(TRUE)
path           <- args[1]
input.file     <- args[2]
filename       <- args[3]
p1.name        <- args[4]
p2.name        <- args[5]
child.name     <- args[6]
perl.to.R.file <- args[7]
graph.output   <- args[8]
verbose        <- args[9]

# Read in tab delimited file of trio, format = SNP Name, chromosome,
# position, parent1, parent1 BAF, parent1 LRR, parent2, parent2 BAF,
# parent2 LRR, child, child BAF, child LRR

# Open data file from perl
#if (verbose) cat("Reading abnormal regions file\n")
abn.regions <- read.delim(file=perl.to.R.file, colClasses=c(rep("numeric", 3),
                          rep("character",2), rep("numeric", 11)))
colnames(abn.regions) <- c("Chr", "Start", "Stop", "Type", "Parent", "Num.SNPs",
                           "Inf.SNPs", "Num.Bases", "Radius", "Midpoint",
                           "Child.Mean.mBAF", "Child.Mean.LRR", "P1.Mean.mBAF",
                           "P1.Mean.LRR", "P2.Mean.mBAF", "P2.Mean.LRR")

#if (verbose) cat("Reading input file\n")
data <- read.delim(file=input.file, header=TRUE, comment.char="",
                   colClasses=c(rep("factor", 2), "numeric", "factor",
                   rep("numeric", 2), "factor", rep("numeric", 2),
                   "factor", rep("numeric", 2)))
data <- na.omit(data)
colnames(data)<-c("SNP Name", "Chr", "Pos", "P1.GType", "P1.BAF", 
                  "P1.LRR", "P2.GType", "P2.BAF", "P2.LRR", "Child.GType", 
                  "Child.BAF", "Child.LRR")
data <- data[grep("[YM-]", data$Chr, invert=T), ]
data$Chr <- as.factor(sub("X", 23, data$Chr))

PlotAbnChr <- function(chromosome)
{
  if (verbose) cat("Plotting Chromosome",chromosome,"\n")
  curr.regions=abn.regions[which(abn.regions$Chr == chromosome),]
  curr.data=data[which(data$Chr == chromosome),]

  moving.avg <- SMA(curr.data$Child.LRR, n=50)  
  
  if (chromosome == 23) chromosome <- "X" 
  x.label <- paste("Chromosome",chromosome, sep="")
  
  if (graph.output == "png" | graph.output == "both") {
    #Create PNG plot
    png(file = paste(path, "/", child.name, "_Chr_", chromosome, 
      ".png", sep=""), width=11, height=8.5, units="in", res=300)    
    par(mfrow=c(3,1), mar=c(0.75, 11, 5, 4),family="sans", font=2, las=0)
    plot(curr.data$Pos, curr.data$Child.LRR,
      xlab     = "", 
      ylab     = "",
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)), 
      ylim     = c(-2,2),
      xaxt     = "n",
      yaxt     = "n",
      pch      = 16, 
      cex      = 0.3,
      font.lab = 2,
      col      = "gray32")
    title(main=child.name, col.main="black", font.main=2, cex.main=4, line=1)
    axis(2,cex.axis=2.5, las=2)
    mtext("LRR", side=2, line=5, font=2, cex=2.5)

    par(new=TRUE)
       
    plot(curr.data$Pos, moving.avg,
      xlab     = "", 
      ylab     = "",
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)), 
      ylim     = c(-2,2), 
      xaxt     = "n",
      yaxt     = "n",
      pch      = 16, 
      cex      = 0.25,
      type     = "l",
      lwd      = 2,
      col      = "green")
    abline(h=0, lty=1, col="black", lwd=2)
    
    par(mar=c(2.875, 11, 2.875, 4))
    
    plot(curr.data$Pos, curr.data$Child.BAF, 
      xlab     = "", 
      ylab     = "",
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)), 
      ylim     = c(0,1),
      xaxt     = "n",
      yaxt     = "n", 
      pch      = 16, 
      cex      = 0.3,
      #cex.lab  = 2.5,
      font.lab = 2,
      col      = "blue")
    axis(2,cex.axis=2.5, las=2)
    mtext("triPOD Results", side=1, line=2.25, font=2, cex=2)
    mtext("BAF", side=2, line=5, font=2, cex=2.5)
  
    par(mar=c(5, 11, 0.75, 4))
      
    plot(1,1, 
      xlab     = "",   
      ylab     = "", 
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)),
      ylim     = c(0,100), 
      xaxt     = "n",    
      yaxt     = "n",
      pch      = 16, 
      cex      = 0.6,
      col      = "white")
    axis(1,cex.axis=2.5)
    mtext("Father", side=2, line=0.5, at=96, font=2, cex=2, las=2)
    mtext("Mother", side=2, line=0.5, at=5, font=2, cex=2, las=2)
    mtext(x.label, side=1, line=3.5, font=2, cex=1.75, las=1)

    for (i in 1:length(curr.regions$Chr)) {
      if (curr.regions$Type[i] == "DEL") color <- "red"
      else if (curr.regions$Type[i] == "AMP") color <- "blue"
      else if (curr.regions$Type[i] == "HD") color <- "gray"
      else if (curr.regions$Type[i] == "UPhD") color <- "purple"
      else if (curr.regions$Type[i] == "UPiD") color <- "green"
      else color <- "black"
    
      min.plotsize = (max(curr.data$Pos) - min(curr.data$Pos)) * 0.001
      if ((curr.regions$Stop[i] - curr.regions$Start[i]) < min.plotsize) {
        radius <- min.plotsize / 2 
      }
      else radius <- curr.regions$Radius[i]
      midpoint <- curr.regions$Midpoint[i]

      if (substr(curr.regions$Parent[i],0,2) == "Mo") {
        roundrect(mid=c(midpoint, 5),radius,4,col=color,lcol=color,lwd=0)
      }
      else if (substr(curr.regions$Parent[i],0,2) == "Fa") {
        roundrect(mid=c(midpoint, 95),radius,4,col=color,lcol=color,lwd=0)
      }
      else {
        roundrect(mid=c(midpoint, 50),radius,4,col=color,lcol=color,lwd=0)
      }
    }
  dev.off() 
  }    
  
  if (graph.output == "pdf" | graph.output == "both") { 
    #Create PDF plot
    pdf(file = paste(path, "/", child.name, "_Chr_", chromosome, ".pdf", sep=""),
      width=11,height=8.5)
    par(mfrow=c(3,1), mar=c(0.75, 11, 5, 4),family="sans", font=2, las=0)

    plot(curr.data$Pos, curr.data$Child.LRR,
      xlab     = "", 
      ylab     = "",
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)), 
      ylim     = c(-2,2),
      xaxt     = "n",
      yaxt     = "n",
      pch      = 16, 
      cex      = 0.3,
      font.lab = 2,
      col      = "gray32")
    title(main=child.name, col.main="black", font.main=2, cex.main=4, line=1)
    axis(2,cex.axis=2.5, las=2)
    mtext("LRR", side=2, line=5, font=2, cex=2.5)

    par(new=TRUE)
       
    plot(curr.data$Pos, moving.avg,
      xlab     = "", 
      ylab     = "",
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)), 
      ylim     = c(-2,2), 
      xaxt     = "n",
      yaxt     = "n",
      pch      = 16, 
      cex      = 0.25,
      type     = "l",
      lwd      = 2,
      col      = "green")
    abline(h=0, lty=1, col="black", lwd=2)
    
    par(mar=c(2.875, 11, 2.875, 4))
    
    plot(curr.data$Pos, curr.data$Child.BAF, 
      xlab     = "", 
      ylab     = "",
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)), 
      ylim     = c(0,1),
      xaxt     = "n",
      yaxt     = "n", 
      pch      = 16, 
      cex      = 0.3,
      #cex.lab  = 2.5,
      font.lab = 2,
      col      = "blue")
    axis(2,cex.axis=2.5, las=2)
    mtext("triPOD Results", side=1, line=2.25, font=2, cex=2)
    mtext("BAF", side=2, line=5, font=2, cex=2.5)
	
    par(mar=c(5, 11, 0.75, 4))
      
    plot(1,1, 
      xlab     = "",   
      ylab     = "", 
      xlim     = c(min(curr.data$Pos), max(curr.data$Pos)),
      ylim     = c(0,100), 
      xaxt     = "n",    
      yaxt     = "n",
      pch      = 16, 
      cex      = 0.6,
      col      = "white")
    axis(1,cex.axis=2.5)
    mtext("Father", side=2, line=0.5, at=96, font=2, cex=2, las=2)
    mtext("Mother", side=2, line=0.5, at=5, font=2, cex=2, las=2)
    mtext(x.label, side=1, line=3.5, font=2, cex=1.75, las=1)

    for (i in 1:length(curr.regions$Chr)) {
      if (curr.regions$Type[i] == "DEL") color <- "red"
      else if (curr.regions$Type[i] == "AMP") color <- "blue"
      else if (curr.regions$Type[i] == "HD") color <- "gray"
      else if (curr.regions$Type[i] == "UPhD") color <- "purple"
      else if (curr.regions$Type[i] == "UPiD") color <- "green"
      else color <- "black"
    
      min.plotsize = (max(curr.data$Pos) - min(curr.data$Pos)) * 0.001
      if ((curr.regions$Stop[i] - curr.regions$Start[i]) < min.plotsize) {
        radius <- min.plotsize / 2 
      }
      else radius <- curr.regions$Radius[i]
      midpoint <- curr.regions$Midpoint[i]

      if (substr(curr.regions$Parent[i],0,2) == "Mo") {
        roundrect(mid=c(midpoint, 5),radius,4,col=color,lcol=color,lwd=0)
      }
      else if (substr(curr.regions$Parent[i],0,2) == "Fa") {
        roundrect(mid=c(midpoint, 95),radius,4,col=color,lcol=color,lwd=0)
      }
      else {
        roundrect(mid=c(midpoint, 50),radius,4,col=color,lcol=color,lwd=0)
      }
    }
  dev.off()  
  } 
}
abn.chr.list <- unique(abn.regions$Chr)
abn.chr.list <- abn.chr.list[which(abn.chr.list <= 23)]
abn.chr.list <- abn.chr.list[order(abn.chr.list)]
tmp<-sapply(abn.chr.list, PlotAbnChr)
cat("The graphics are complete.\n")
q(save = "no", status = 0, runLast = TRUE)
