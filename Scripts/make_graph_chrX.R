
args = commandArgs(trailingOnly=TRUE)

if (length(args)<3) {
  stop("At least 3 arguments must be supplied (patient's name, input file, output file).n", call.=FALSE)
} else {
	patient=args[1]
	file=args[2]
	output=args[3]
	size=args[4]
}

# chromosome names and positoins
chr_name=paste("chr",seq(1,22,1),sep="") 
chr_name=c(chr_name,"chrX")
chr_position=c(124478211,370053186,590297730,784552788,970429194,1146601313,1311677289,1463919594,1605686270,1741782340,1876224362,2010405327,2134225146,2244929169,2349446622,2445611389,2532409282,2614224645,2683720096,2745250987,2800828062,2849592288, 2953021970)
chr_limits=c(0,248956422, 491149951, 689445510, 879660065, 1061198324, 1232004303, 1391350276, 1536488912, 1674883629, 1808681051, 1943767673, 2077042982, 2191407310, 2298451028, 2400442217, 2490780562, 2574038003, 2654411288, 2713028904, 2777473071, 2824183054, 2875001522, 3031042417)

####
pdf(output, width=10, height=3)


plot(c(0,max(chr_limits)),c(0,100),xlim=c(1,max(chr_limits)), yaxs = "i", xaxs = "i",type="n",xlab="", ylab="",axes=F, main=paste("Homozygous Regions for ", patient, "\n Total = ",size," Mb (autosomes)", sep=""))

box()
## For the x axis
axis(1,tck=0, at=chr_position,  labels=chr_name,las=2,cex=0.8)
axis(side = 1, tck = -.10, labels = NA, at=chr_limits)
## For the y axis
#axis(side=2,tck=-.010, at=c(0, 25, 50, 75, 100), labels=c(0, 25, 50, 75, 100),las=1)


## PLOT THE CHROMOSOME LIMIT"S
abline(v=chr_limits[2:23],col='grey',lty=2)


data<-read.table(file,sep="\t",comment.char="@",fill=T) #import the homozygous region file
if(dim(data)[1]>4){
	data<-read.table(file,header=F,sep="\t") #import the homozygous region file
	### PRINT THE HOMOZYGOSITY DATA
	chr_pos<-as.numeric(substr(data$V1,4,5)) #extract the chromosome number
	chr_pos[which(is.na(chr_pos))]=23
	xleft<-data$V2+chr_limits[chr_pos] # extract the left position (position on the chromosome + chromosome limit)
	xright<-data$V3+chr_limits[chr_pos] # extract the right position (position on the chromosome + chromosome limit)
	rect(xleft, 0, xright, 100, col = "blue", border = NA) #print it
}

dev.off()


