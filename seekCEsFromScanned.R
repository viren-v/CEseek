seekCEsFromScanned<-function(CEscan.out, test.names, bg.names, p.cut,  out.dir, output.html="no"){
	
	
	pan.outes <- lapply(CEscan.out[[2]], function(x) x[[1]])
	pan.outed <- lapply(pan.outes, function(x) lapply(x, function(y){ if (is.list(y)){output<-Reduce(c, y[2:length(y)])} else{output=y}; return(output)} ))
	comb.vector <- names(CEscan.out[[2]])
	j.vector <- names(CEscan.out[[2]][[1]][[1]])
	
###############################################################################################	
makePWMpng<- function(x, out.dir, filename){
getProbMatrixfromDNAStringset <- function(x) {
ce.ce.seq.unlist <- x
vkc.data <- t(sapply(ce.ce.seq.unlist, function(x) strsplit(as.character(x), split="")[[1]]))
vkc.data.2 <- apply(vkc.data, 2, function(x) c(length(which(x=="A")), length(which(x=="C")), length(which(x=="G")), length(which(x=="T")))/length(x))
vkc.pfm <- apply(vkc.data, 2, function(x) c(length(which(x=="A")), length(which(x=="C")), length(which(x=="G")), length(which(x=="T"))))
rownames(vkc.pfm)<- c("A", "C", "G", "T")
pbm<- apply(vkc.pfm,2, function(x)x/sum(x))
return(pbm)
}
	filename.new <- gsub("[/]", "_", gsub(":", "_", filename)) #gsub("[:|/]", "_", filename)
	vkc.pbm<- getProbMatrixfromDNAStringset(x)
	out.file.name <- paste(out.dir, filename.new, ".png", sep="")
	png(out.file.name, height=300, width=100+(30*ncol(vkc.pbm)))
	seqLogo(vkc.pbm, xfontsize=10, yfontsize=10)
	dev.off()
	return(vkc.pbm)
}
#################################################################################################	
getFishersPvalues_loop2 <- function(xm, ym, x.le, y.le){
		test.matrix <- cbind(xm, ym, check.names=T)[cbind(xm, ym, check.names=T)[,3]==1,1:2]
	rm(xm, ym)
	depletion.enrich.x.fisher  <- matrix(,nrow(test.matrix), 2)
	for(i in 1:nrow(test.matrix)){
		x <- test.matrix[i,]
	depletion.enrich.x.fisher[i,] <- c(fisher.test(matrix(c(x[1], x[2], x.le-x[1], y.le-x[2] ),2,2) , alternative="greater")$p.value, fisher.test(matrix(c(x[1], x[2], x.le-x[1], y.le-x[2] ),2,2) , alternative="less")$p.value)
	
	}
	colnames(depletion.enrich.x.fisher) <- c("p.enrichment", "p.depletion")
	percent.seqs.fold <- cbind(test.matrix, (test.matrix[,1]*100)/x.le, (test.matrix[,2]*100)/y.le, log2(((test.matrix[,1])/x.le)/((test.matrix[,2])/y.le)))
	
	colnames(percent.seqs.fold) <- c("# of Test seqs", "# of background seqs", "% of Seq with CE in vkc.sample", "% of Seq with CE in Background", "log2 Enrichment")
	
	output.fish <- data.frame(stats=cbind(depletion.enrich.x.fisher, percent.seqs.fold) )
		
	return(output.fish)
	}

################################################################################################
getFishersPvalues_loop <- function(xm, ym, x.le, y.le){
	## xm and ym are zoops matrix
	# first convert it into a single column each
	c1 <- rep(colnames(xm), each=nrow(xm))
	c2 <- rep(rownames(xm), times=ncol(xm))
	test.matrix <- cbind(c(xm), c(ym))
	rm(xm, ym)
	depletion.enrich.x.fisher  <- matrix(,nrow(test.matrix), 2)
	for(i in 1:nrow(test.matrix)){
		x <- test.matrix[i,]
	depletion.enrich.x.fisher[i,] <- c(fisher.test(matrix(c(x[1], x[2], x.le-x[1], y.le-x[2] ),2,2) , alternative="greater")$p.value, fisher.test(matrix(c(x[1], x[2], x.le-x[1], y.le-x[2] ),2,2) , alternative="less")$p.value)
	#print(i)
	}
	colnames(depletion.enrich.x.fisher) <- c("p.enrichment", "p.depletion")
	percent.seqs.fold <- cbind(test.matrix, (test.matrix[,1]*100)/x.le, (test.matrix[,2]*100)/y.le, log2(((test.matrix[,1])/x.le)/((test.matrix[,2])/y.le)))
	
	colnames(percent.seqs.fold) <- c("# of Test seqs", "# of background seqs", "% of Seq with CE in vkc.sample", "% of Seq with CE in Background", "log2 Enrichment")
	
	output <- data.frame(TFs=c1, structure=c2, stats=cbind(depletion.enrich.x.fisher,percent.seqs.fold) )
		
	return(output)
	
}



################################################################################################

	test.le <- length(test.names)
	bg.le <- length(bg.names)
library(seqLogo)	
pan.out <- lapply(pan.outed, function(x) x[intersect(names(x), j.vector)])
t.ct <- lapply(pan.out, function(x) lapply(x, function(y){ test.names.x <- intersect(names(y), test.names); return(y[test.names.x]) }))
b.ct <- lapply(pan.out, function(x) lapply(x, function(y){ bg.names.x <-   intersect(names(y), bg.names); return(y[bg.names.x])    }))	

t.ctsum <- sapply(t.ct, function(x) length(Reduce(union, sapply(x, names))))
b.ctsum <- sapply(b.ct, function(x) length(Reduce(union, sapply(x, names))))
out.stats.2 <- getFishersPvalues_loop2(t.ctsum, b.ctsum, test.le, bg.le)




#print("Stage1 complete")
tce.list <- lapply(t.ct, function(x) sapply(x, function(y) length(unique(names(y)))))
bce.list <- lapply(b.ct, function(x) sapply(x, function(y) length(unique(names(y)))))
print("Calculating significance")
tce <- matrix(0,length(j.vector), length(comb.vector))
rownames(tce)<- j.vector
colnames(tce)<- comb.vector
for(i in 1:length(comb.vector)){tce[names(tce.list[[i]]),i]<- tce.list[[i]]}
#print("3")
bce <- matrix(0,length(j.vector), length(comb.vector))
rownames(bce)<- j.vector
colnames(bce)<- comb.vector
for(i in 1:length(comb.vector)){bce[names(bce.list[[i]]),i]<- bce.list[[i]]}
#print(head(tce))
out.stats <- getFishersPvalues_loop(tce, bce, test.le, bg.le)
#print(head(out.stats))
#out.stats <- out.stats[,c(1:4,9,5:8)]

rownames(out.stats)<- paste(out.stats[,1], out.stats[,2], sep="_")
#print("4")
t.ct.nz <- lapply(t.ct, function(x) x[sapply(x, length)>0])
b.ct.nz <- lapply(b.ct, function(x) x[sapply(x, length)>0])

t.ct.unlist.nz <- unlist(t.ct.nz)
b.ct.unlist.nz <- unlist(b.ct.nz)
names(t.ct.unlist.nz) <- paste(rep(names(t.ct.nz), times=sapply(t.ct.nz, length)), unlist(sapply(t.ct.nz, names)), sep="_")
names(b.ct.unlist.nz) <- paste(rep(names(b.ct.nz), times=sapply(b.ct.nz, length)), unlist(sapply(b.ct.nz, names)), sep="_")

stat.names<- paste(out.stats[,1], out.stats[,2], sep="_")


t.ct.unlist.nz.signi <- t.ct.unlist.nz[rownames(out.stats)[out.stats[,3]<=p.cut]]
b.ct.unlist.nz.signi <- b.ct.unlist.nz[rownames(out.stats)[out.stats[,4]<=p.cut]]

print(paste("Count of CEs enriched in test set= ", length(t.ct.unlist.nz.signi), sep=""))
print(paste("Count of CEs depleted in test set= ", length(b.ct.unlist.nz.signi), sep=""))

dir.create(out.dir)
dir.create(paste(out.dir,"/testLogo/", sep=""))
dir.create(paste(out.dir,"/backgroundLogo/", sep=""))
out.dir.test <- paste(out.dir,"/testLogo/", sep="")
out.dir.bg <- paste(out.dir,"/backgroundLogo/", sep="")

#print("5")
t.ct.unlist.nz.signi.prob <- list()
b.ct.unlist.nz.signi.prob <- list()

if(length(t.ct.unlist.nz.signi)>0){
for(i in 1:length(t.ct.unlist.nz.signi)){t.ct.unlist.nz.signi.prob[[i]] <- makePWMpng(t.ct.unlist.nz.signi[[i]], out.dir.test, names(t.ct.unlist.nz.signi)[i])}
names(t.ct.unlist.nz.signi.prob) <- names(t.ct.unlist.nz.signi)
} else {}
if(length(b.ct.unlist.nz.signi)>0){
for(i in 1:length(b.ct.unlist.nz.signi)){b.ct.unlist.nz.signi.prob[[i]] <- makePWMpng(b.ct.unlist.nz.signi[[i]], out.dir.bg, names(b.ct.unlist.nz.signi)[i])}
names(b.ct.unlist.nz.signi.prob) <- names(b.ct.unlist.nz.signi)
} else{}
stat.names.signi <- c(names(t.ct.unlist.nz.signi), names(b.ct.unlist.nz.signi))

stat.names.signi.forBrowser <- gsub("[/]", "_", gsub(":", "_", stat.names.signi))
print(stat.names.signi.forBrowser)
#print(head(stat.names.signi.forBrowser))
test.seq.logo.path = paste('<img src="', "./testLogo/", stat.names.signi.forBrowser, '.png">', sep='')
bg.seq.logo.path = paste('<img src="', "./backgroundLogo/", stat.names.signi.forBrowser, '.png">', sep='')

stats.final <- out.stats[stat.names.signi,]
output<-cbind(testLogo=test.seq.logo.path, backgroundLogo=bg.seq.logo.path, stats.final)
if(output.html=="yes"){
library(SortableHTMLTables)

sortable.html.table(output, "CE.results.html", output.directory=out.dir, page.title = out.dir)
} else {}
	
final.output<- list(TestProbMatrix=t.ct.unlist.nz.signi.prob, BackgroundProbMatrix=b.ct.unlist.nz.signi.prob, CEEnrichment.Stats=output, CE.All.Stats= out.stats,TFcombn.Enrichment=out.stats.2)

return(final.output)		
}