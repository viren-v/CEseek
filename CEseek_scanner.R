CEseek_scanner <- function(input.seqs, tf.pwms, length.gap, half.gap, match.thres) {
le.gap <- length.gap
#########################################################################################
iranges2df<-function(my.object, char){	
	a <- as.data.frame(as(my.object, "IRanges"))
	b <- as.data.frame(as.character(my.object))
	c <- matrix(char, nrow(a), 1)
	output<- cbind(a,b,c)
	colnames(output)<- c("start", "end", "width", "Sequence", "Strand Relative to TSS")
	return(output)
}

################################################################################################

scan_TFpwms_list_on_a_sequence <- function(matrices, promoter.current, match.thres){
## match.thres is a list of PWM associated stringency and has the same structure as the pwms have
curTFpwms.on.promoter.currents <- list()
vkc.names <- character()
jango <- sapply(matrices, length)
for (i in 1:length(matrices)){
curTFpwm.set <- matrices[[i]]	
match.thres.set <- match.thres#[[i]]
for (j in 1:jango[i]){
	
curTFpwm <- curTFpwm.set[[j]]
cur.match.thres <- match.thres.set#[[j]]
curTFpwm.on.promoter.current.forward <- matchPWM(curTFpwm, promoter.current[[1]], as.character(cur.match.thres[[1]]))
curTFpwm.on.promoter.current.reverse <- matchPWM(reverseComplement(curTFpwm), promoter.current[[1]], as.character(cur.match.thres[[1]]))

forward.score.raw <- PWMscoreStartingAt(curTFpwm, promoter.current[[1]],starting.at = start(curTFpwm.on.promoter.current.forward))
reverse.score.raw <- PWMscoreStartingAt(reverseComplement(curTFpwm), promoter.current[[1]], starting.at = start(curTFpwm.on.promoter.current.reverse))
forward.score <- forward.score.raw/maxScore(curTFpwm)
reverse.score <- reverse.score.raw/maxScore(reverseComplement(curTFpwm))

curTFpwm.on.promoter.current.forward.df <- cbind(iranges2df(curTFpwm.on.promoter.current.forward, "+"), forward.score)
curTFpwm.on.promoter.current.reverse.df <- cbind(iranges2df(curTFpwm.on.promoter.current.reverse, "-"), reverse.score)

colnames(curTFpwm.on.promoter.current.forward.df)<-c("start", "end", "width", "Sequence", "Strand Relative to TSS", "Score")
colnames(curTFpwm.on.promoter.current.reverse.df)<-c("start", "end", "width", "Sequence", "Strand Relative to TSS", "Score")

curTFpwm.on.promoter.current<- rbind(curTFpwm.on.promoter.current.forward.df, curTFpwm.on.promoter.current.reverse.df)
## colnames(curTFpwm.on.promoter.current) <- c("start", "end", "width", "sequence", "Strand Relative to TSS", "Score")

curTFpwms.on.promoter.currents[[length(curTFpwms.on.promoter.currents)+1]]<- curTFpwm.on.promoter.current

vkc.name <- paste(names(matrices)[i], j)
vkc.names[length(vkc.names) + 1] <- vkc.name 

}
}
names(curTFpwms.on.promoter.currents)<-vkc.names
return(curTFpwms.on.promoter.currents)
}
################################################################################################

scan_TFpwms_list_on_sequence_list <- function(mdb.list, seq.list, match.thres){
mdb.nonzero <- mdb.list[(!sapply(mdb.list, length)==0)]

curTFpwms.on.promoters.current <- list(list())
# d.tss<- 100000
for (i in 1:length(seq.list)){
seq.list.curr <- seq.list[i]
curTFpwms.on.promoter.current <- scan_TFpwms_list_on_a_sequence(mdb.nonzero, seq.list.curr, match.thres)

curTFpwms.on.promoters.current[[i]]<- curTFpwms.on.promoter.current
}
names(curTFpwms.on.promoters.current)<-names(seq.list)
return(curTFpwms.on.promoters.current)
}

################################################################################################
convert_into_new_nested_list_format<- function(scan_out.cur, pwms){
	test <- names(scan_out.cur)
	test2 <- sub("[ ]+[0-9]*", "", test)
	pwm_groups <- match(test2, unique(test2))
	mynames <- unique(test2)
	my.list <- list(); 
	for (i in 1:length(mynames)){
		pe <- which(pwm_groups==i);
		my.list[[i]]<- list()
		for(j in 1:length(pe)){
			my.list[[i]][[j]]<- scan_out.cur[[pe[j]]]
	}
	
	}
	
	names(my.list)<- names(pwms)
	for(i in 1:length(my.list)){ names(my.list [[i]] ) <- names(pwms[[i]]) }
	return(my.list)
}


################################################################################################

olRanges.igStrand.wOrientRange_dbud<- function(queryO, subjectO, halfgap, scanned.Seq, seq.frame.list, seq.name) {
       TFs1 <- names(queryO) 
       TFs2 <- names(subjectO)
            queryO <- queryO[[1]]
            subjectO <- subjectO[[1]]  
            queryO<- GRanges("scanned.Seq", IRanges(queryO[,1], queryO[,2]), strand=queryO[,5])
            subjectO<- GRanges("scanned.Seq", IRanges(subjectO[,1], subjectO[,2]), strand=subjectO[,5])
        query <- promoters(queryO, upstream=halfgap, downstream= halfgap+ width(queryO)[1])
        subject <- promoters(subjectO, upstream=halfgap, downstream= halfgap+ width(subjectO)[1])
        
        polix <- cbind(c("+","+","+","+","-","-","-","-"), c("+","+","-","-","+","+","-","-"), c("+","-","+","-","+","-","+","-"))
        poliz <- paste(polix[,1], polix[,2], polix[,3], sep="|")
        mrig <- cbind(c("+","-","+","-","-","+","-","+"), c("+","+","-","-","-","-","+","+"))
        mrigz <- paste(mrig[,1], mrig[,2], sep="|")
       	names(mrigz)<- poliz
        
	olindex <- as.matrix(findOverlaps(query, subject, maxgap = 1L, ignore.strand=TRUE))
	if(nrow(olindex)>0){
	ad.sign <- array(, nrow(olindex))
	ad.sign[start(subjectO[olindex[,2]])-start(queryO[olindex[,1]])<0]<- "+"
	ad.sign[!start(subjectO[olindex[,2]])-start(queryO[olindex[,1]])<0]<- "-"
	
	nolindex <- cbind(olindex, ad.sign)
	}
	
	if(nrow(olindex)>0){
        query <- query[olindex[,1]]
        subject <- subject[olindex[,2]]
        query.o <- queryO[olindex[,1]]
        subject.o <- subjectO[olindex[,2]]
                
        olma <- cbind(Qstart=start(query), Qend=end(query), Sstart=start(subject), Send=end(subject))
        iolma<- cbind(Qstart=start(query.o), Qend=end(query.o), Sstart=start(subject.o), Send=end(subject.o))
		Q.ir <- IRanges(start(query), end(query))
		S.ir <- IRanges(start(subject), end(subject))
		dis.o <- distance(query.o, subject.o, ignore.strand=TRUE)
	#	dis.o[dis.o==0]<- -(width(setdiff(query.o[dis.o==0], subject.o[dis.o==0])))
		if(length(which(dis.o==0))==0){
			configs.SJ <- cbind(as.matrix(as.data.frame(query.o)$strand), as.matrix(as.data.frame(subject.o)$strand), as.matrix(nolindex)[,3])
			 }
		else if( length(which(dis.o==0)) ==1){
			x<-iolma[dis.o==0,]
		dis.o[dis.o==0]<- -length(intersect(x[1]:x[2], x[3]:x[4]))
		configs.SJ <- cbind(as.matrix(as.data.frame(query.o)$strand), as.matrix(as.data.frame(subject.o)$strand), as.matrix(nolindex)[,3])
			
		
			
			}	else {
		dis.o[dis.o==0]<- apply(iolma[dis.o==0,], 1, function(x) -length(intersect(x[1]:x[2], x[3]:x[4])))
		configs.SJ <- cbind(as.matrix(as.data.frame(query.o)$strand), as.matrix(as.data.frame(subject.o)$strand), as.matrix(nolindex)[,3])
				
		}
		configs.p <- paste(configs.SJ[,1],configs.SJ[,2], configs.SJ[,3], sep="|")
			configs.SJP <- mrigz[configs.p]
			
			gr.frSeq <- t(apply(iolma, 1, range))
			rownames(gr.frSeq) <- configs.p
			#rcd.n <- 
			#rcd.z <- 
			
			
		output<- cbind(dis.o, configs.SJP, gr.frSeq)
		
		
		rownames(output) <- paste(output[,1], output[,2], sep="|")
		my.seq.list<- list() 
		for(i in 1:nrow(gr.frSeq)){
			if (length(intersect(rownames(gr.frSeq)[i], names(mrigz[5:8])))==0) { my.seq.list[[i]] <- DNAStringSet(subseq(scanned.Seq, start=min(gr.frSeq[i,]), end=max(gr.frSeq[i,]))); names(my.seq.list[[i]])<- rep(seq.name, length(my.seq.list[[i]])) } else {my.seq.list[[i]] <- DNAStringSet(reverseComplement(subseq(scanned.Seq, start=min(gr.frSeq[i,]), end=max(gr.frSeq[i,])))); names(my.seq.list[[i]])<- rep(seq.name, length(my.seq.list[[i]])) }
		
		}
		rownames(gr.frSeq) <- paste(output[,1], output[,2], sep="|")
		names(my.seq.list) <- rownames(gr.frSeq)
		 
		for(i in 1:length(seq.frame.list)){
		seq.frame.list[[i]] <- my.seq.list[which(names(my.seq.list)==names(seq.frame.list)[i])]
		}
		if(sum(sapply(seq.frame.list, length))==0) {final.output<- list(scanned.Seq)} else{ final.output<- list(output, seq.frame.list, scanned.Seq)}
		#return(final.output)
		}	else {
			final.output<- list(scanned.Seq)
			}
		#print("hi")
		return(final.output)}



################################################################################################

makeCEframe <- function(scaned_out, halfgap, pwms, scanned.seqs, le.gap){

	scan_out <- lapply(scaned_out, function(x) convert_into_new_nested_list_format(x, pwms))
	
	pos.c <- c("+|+", "-|-", "+|-", "-|+")
	
	my.frame.list<- lapply(as.list(le.gap[1]:le.gap[2]), function(x) paste(x, pos.c, sep="|"))
	my.frame <- unlist(lapply(as.list(le.gap[1]:le.gap[2]), function(x) paste(x, pos.c, sep="|")))
	large.frame <- rep(0, length(my.frame))
	large.frame.list <- vector(mode="list", length=length(my.frame))
	names(large.frame)<- my.frame
	names(large.frame.list)<- my.frame
	
	tfs <- names(pwms)
	c2 <- combn(tfs, 2)
	 c2.list <- lapply(apply(c2,2, as.list), unlist)
	pwms.names<- sapply(pwms, names)
	pwms.alt<- sapply(pwms, names)
	
	for(i in 1:length(pwms.alt)){for(j in 1:length(pwms.alt[[i]])){pwms.alt[[i]][[j]]<- paste(names(pwms.alt)[i], j, sep="_")}}
	# names(pwms.names)<- unlist(pwms.alt)
	
	c2.2 <- apply(c2, 2, function(x) {output<-expand.grid(pwms.names[x]); if(ncol(output)==1){output<- matrix(t(output),1,2)} else {output=output}; return(output)})
	if(class(c2.2)=="matrix"){c2.2 <- lapply(apply(c2.2,2, list), function(x)  matrix(unlist(x),1,2))} else {}
	
	
	c2.2.frame <- c2.2
 	for(i in 1:length(c2.2.frame)){le<- nrow(c2.2.frame[[i]]); c2.2.frame[[i]]<- list(); for(j in 1:le){c2.2.frame[[i]] [[j]]<- as.list(large.frame)}; names(c2.2.frame[[i]])<- paste(c2.2[[i]][,1], c2.2[[i]][,2], sep=".w." )}
 	
 	c2.2.frame.seqs <- c2.2
 	for(i in 1:length(c2.2.frame.seqs)){le<- nrow(c2.2.frame.seqs[[i]]); c2.2.frame.seqs[[i]]<- list(); for(j in 1:le){c2.2.frame.seqs[[i]] [[j]]<- as.list(large.frame)}; names(c2.2.frame.seqs[[i]])<- paste(c2.2[[i]][,1], c2.2[[i]][,2], sep=".w." )}

 
	names(c2.2.frame) <- paste(c2[1,], c2[2, ], sep=".W.")
	names(c2.2.frame.seqs) <- paste(c2[1,], c2[2, ], sep=".W.")
	
	
	
	
	for(i in 1:length(scan_out)){
		
		scan_out.cur<- scan_out[[i]]
		scanned.seqs.cur <-scanned.seqs[[i]] 
		scanned.seqs.cur.name <- names(scanned.seqs)[i]
		
		#print(i)
				scan.pwms.m <- pwms.names[which(sapply(lapply(scan_out.cur, function(x) which(sapply(x, function(y)nrow(y)>0))), length)>0)]
		scan.pwms <- names(scan.pwms.m)[!is.na(names(scan.pwms.m))]
		scan.pwms.c2.2.id <- which(sapply(c2.list, function(x) length(intersect(x, scan.pwms)))==2)
		
		if(length(scan.pwms.c2.2.id)>0)	{
		for(k in 1:length(scan.pwms.c2.2.id)){
			TF1 <- c2.list[[scan.pwms.c2.2.id[k] ]][1]
			TF2 <- c2.list[[scan.pwms.c2.2.id[k] ]][2]
			
			scan.tf12 <- scan_out.cur[c(TF1, TF2)]
			cur.3.id <- sapply(scan.tf12, function(x) names(x)[which(sapply(x, nrow)>0)])
			if(length(unique(cur.3.id))>1){
				
						ce.scans.result <- olRanges.igStrand.wOrientRange_dbud(scan_out.cur[[TF1]], scan_out.cur[[TF2]],  halfgap, scanned.seqs.cur, large.frame.list, scanned.seqs.cur.name )	
			if(length(ce.scans.result)>1){
			
			ce.scans.results <- ce.scans.result[[1]]
			ce.scans.results.collapsed <- lapply(as.list(table(rownames(ce.scans.results))), function(x) rep(scanned.seqs.cur.name, as.numeric(x)))
			
			
			
			binned.info <- ce.scans.results.collapsed #
			retrieved.seqs <- ce.scans.result[[2]]
			binned.seqs <- retrieved.seqs[sapply(retrieved.seqs, length)>0]

			
			for(clu in 1:length(binned.info)) {
			
			
			c2.2.frame.seqs[[scan.pwms.c2.2.id[k] ]][[1]][[names(binned.seqs)[clu]]] <- c(c2.2.frame.seqs[[scan.pwms.c2.2.id[k] ]][[1]][[names(binned.seqs)[clu]]], binned.seqs[[names(binned.seqs)[clu]]])
			
			}
			}
		}	
		}
		
		}
			
#print(i); print(Sys.time()) 
}


output.finalized<- 	list(c2.2.frame, c2.2.frame.seqs)

return(output.finalized)

}
out.1 <- scan_TFpwms_list_on_sequence_list(tf.pwms, input.seqs, match.thres)
out.2 <- makeCEframe(out.1, half.gap, tf.pwms, input.seqs, le.gap)

return(out.2)
}