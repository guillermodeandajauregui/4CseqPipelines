#########################################
#Splinter et al. 2012 Methods
#Supplemental R script
#
#Within code is contained that can generate spider plots.

#dummy example
#> source("spider_plots.R")
#> x1 <- seq(0,110e6, by=10e6)
#> x2 <- seq(0,110e6, by=10e6)+2e6
#> makeSpiderGramSingle(dom=cbind(x1,x2), chrom.len=150e6, vp.loc=130e6, col=rgb(0,0,1)) 

#The genes you have to add your self since these are different for every species and chromosome.

#written by Elzo de Wit (2012)
#mail: e.wit@hubrecht.eu
#Hubrecht Institute
#Utrecht, the Netherlands
#########################################


#draw a tubular shaped chromosome
drawChrom <- function(chrom.len = 50e6, max.wid = 100e6, hei=10, chrom.wid=1, y.loc=5){
	r <- seq(0,2*pi, len=1000)
	chrom.wid = chrom.wid/2

	dim.val <- par("din")
	left  <- r[501:1000]
	right <- r[1:500]
	
	correction <- ((max.wid*chrom.wid)/(2*hei)) / (dim.val[1]/dim.val[2])

	x <- sin(left)*correction + correction
	y <- cos(left)*chrom.wid+y.loc

	x1 <- x
	y1 <- y

	x <- c(x, sin(right)*correction+chrom.len -correction)
	y <- c(y, cos(right)*chrom.wid+y.loc )

	x2 <- sin(right)*correction+chrom.len -correction
	y2 <- cos(right)*chrom.wid+y.loc
	polygon(x,y, col='white')
	
	col.seq <- seq(0.5,1,len=250)
	col.seq <- c(col.seq,rev(col.seq))
	cols <- rgb(col.seq,col.seq,col.seq)
	segments(x1,y1,x2,rev(y2), col = cols)
	polygon(x,y, lwd=2)

}	


#function for drawing two chromosomes
drawLocalChrom <- function(labels=c("Xa","Xi"), yloc1 = 3, yloc2 = 5, wid = 2, num.chrom=1, chrom.len){
	plot(c(0,chrom.len), c(yloc1-wid,yloc2+wid), type='n', axes=F, xlab="chromosomal location", ylab="", cex.lab=1)

	mb <- floor(chrom.len/1e6)
	mb10 <- floor(mb/10)*10
	#draw an axis
	label <- c(seq(0,mb10, by=10),mb)
	segments(label*1e6, -1e9, label*1e6, 1e9, lwd=2, lty=2, col='grey90')
	mid.y <- (yloc1+yloc2)/2
	segments(label*1e6, mid.y-0.15, label*1e6, mid.y+0.15, lwd=2)
	segments(0, mid.y, mb*1e6, mid.y, lwd=2)


	#active X
	if(num.chrom==2){
		drawChrom(chrom.len=chrom.len, max.wid=chrom.len, y.loc=yloc1)
		text(-1e6,yloc1, labels[2],cex=1)
	}	
	#inactive X 
	drawChrom(chrom.len=chrom.len, max.wid=chrom.len, y.loc=yloc2)
	
	text(-1e6,yloc2, labels[1],cex=1)
	at <- 0:mb
	axis(1, at=at*1e6, labels=NA, cex.axis=1, lwd=1,las=1)
	axis(1, at=label*1e6, labels=label, cex.axis=1, lwd=2,las=1)
	axis(3, at=at*1e6, labels=NA, cex.axis=1, lwd=1,las=1)
	axis(3, at=label*1e6, labels=label, cex.axis=1, lwd=2)
}

#draw the splines showing the interactions
drawSplines.domain <- function( dom, vp.loc, y.base, y.arc, plot=F, col='black', chrom.size=166e6, relative=F){
	if(nrow(dom) == 0)
		return
	for(i in 1:nrow(dom)){
		#start <- min(dom[i,1], dom[i,2])
		#end <- max(dom[i,1], dom[i,2])
		start <- dom[i,1]
		end <- dom[i,2]
		xspline(c(vp.loc, (end + vp.loc)/2, end, start, (end + vp.loc)/2, vp.loc), c(y.base,y.arc,y.base,y.base,y.arc,y.base), open=F, shape=c(0,1,0,0,1,0), col=col, border=col)
	}
}

#dom:       matrix or data.frame with two columns, start and end position of the interactions
#vp.loc:    the position of the viewpoint
#color:     color of the interaction splines
#gene:      data.frame containing the columns for the start, end and strand of gene
#labels:    what should be put on the left side of the plot
#chrom.len: the length of the chromosome

makeSpiderGramSingle <- function( dom, vp.loc, color='black', gene="", labels=c(""), chrom.len ){

	#draw two chromosomes
	drawLocalChrom( labels = labels, num.chrom=1, chrom.len = chrom.len )

	#and the splines
	drawSplines.domain(dom, vp.loc=vp.loc, plot=F, relative=F, y.arc=8, y.base=5.5, col=color)
		
	#draw the genes
	if(! is.null(nrow(gene))){
		gene <- gene[gene[,2]-gene[,1] < 2e5,]
		rect(gene[,1],5, gene[,2], ifelse(gene[,3]=='+',5.5,4.5), col='black')
	}	

}

