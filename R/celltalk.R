# mycolors <- hue_pal(c=100)(25)

#' To present a circos plot
#' @param Interact Interact list returned by findLRpairs
#' @param ident To highlight the interaction of a specific identity class between others; if 'NULL', plot interaction for all identity classes
#' @return Circos plot showing the ligand-receptor interaction
#' @export
circosPlot <- function(Interact,ident=NULL){
	options(stringsAsFactors=F)
	Interact.num.dat <- Interact$InteractNumer
	# Interact.num.dat = Interact.num.dat[Interact.num.dat$LR.Number!=0,]
	# cluster <- unique(c(Interact.num.dat$Cell.From,Interact.num.dat$Cell.To))
	# mycolors <- hue_pal(c=100)(length(cluster))
	circos.clear()
	circos.par(start.degree=90, clock.wise=F)
	if (is.null(ident)){
		chordDiagram(Interact.num.dat,annotationTrack = c("name","grid"),transparency=0.1)

		}else{

			cols <- rep('gray',nrow(Interact.num.dat))
			cols[(Interact.num.dat$Cell.From==ident) | (Interact.num.dat$Cell.To==ident)] <- 'red'
			chordDiagram(Interact.num.dat,annotationTrack = c("name","grid"),transparency=0.1,col=cols)

		}
}

#' To find marker ligands and marker receptors
#' @param expr.mat Matrix or data frame of expression matrix, with genes in rows and cells in columns
#' @param lable  Vector of identity classes of cells in the expression matrix
#' @param method Method used for differential expression test, either 'wilcox.test' or 't.test'
#' @param p.adjust Method used for p value correction for multiple differential expression test; see p.adjust function for more information
#' @return Data frame containing the differential expression test
#' @export
findDEGs <- function(expr.mat, lable, method='wilcox.test', p.adjust='BH'){
	if (!is.factor(lable)){ lable <- as.factor(lable) }
	ident.level <- levels(lable)
	
	lable <- as.character(lable)
	### p.value calculated
	if (method=='wilcox.test'){
		test.all.res <- lapply(ident.level, function(each.level) {
			print(paste0('Identifying marker genes for cluster ',each.level,' ...'))
			cell.ident <- which(lable==each.level)
			cell.other <- which(lable!=each.level)

			test.res <- apply(expr.mat, 1, function(row.expr){
				logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
				p.value <- wilcox.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value
				c(logFC, p.value)
			})

			test.res <- as.data.frame(t(test.res))
			colnames(test.res) <- c('avg_log2FC','p_val')
			test.res$p_val_adj <- p.adjust(test.res$p_val, method=p.adjust)
			test.res$cluster <- each.level
			test.res$gene <- rownames(test.res)
			return(test.res)

		})
		test.all.res <- do.call(rbind, test.all.res)
	}else if(method=='t.test'){
		test.all.res <- lapply(ident.level, function(each.level) {
			print(paste0('Identify marker genes for ',each.level,' ...'))
			cell.ident <- which(lable==each.level)
			cell.other <- which(lable!=each.level)

			test.res <- apply(expr.mat, 1, function(row.expr){
				logFC <- log(mean(expm1(row.expr[cell.ident])) +1, base=2) - log(mean(expm1(row.expr[cell.other])) +1, base=2)
				p.value <- t.test(x=row.expr[cell.ident], y=row.expr[cell.other])$p.value
				c(logFC, p.value)
			})

			test.res <- as.data.frame(t(test.res))
			colnames(test.res) <- c('avg_log2FC','p_val')
			test.res$p_val_adj <- p.adjust(test.res$p_val, method=p.adjust)
			test.res$cluster <- each.level
			test.res$gene <- rownames(test.res)
			return(test.res)

		})
		test.all.res <- do.call(rbind, test.all.res)
	}else{
		stop("select t.test or wilcox.test to conduct differential analysis")
	}

	return(test.all.res)
}





#' To find marker ligands and marker receptors
#' @param marker.dat Data frame containing information of marker genes
#' @param species species, either 'hsapiens', 'mmusculus', or 'rnorvegicus' 
#' @param logFC.thre logFC threshold, marker genes with a logFC > logFC.thre will be considered
#' @param p.thre p threshold, marker genes with a adjust p value < p.thre will be considered
#' @return List containing the ligand-receptor interaction information
#' @export
findLRpairs <- function(marker.dat, species, logFC.thre=0.25, p.thre=0.01){
	options(stringsAsFactors=F)
	data('cellTalkData', package=('celltalk'))

	lr.pair.dat <- cellTalkData$DataLR[[species]]
	ligs <- lr.pair.dat$L
	reps <- lr.pair.dat$R

	cluster.level <- as.factor(unique(marker.dat$cluster))
	num.cluster <- length(cluster.level)
	Interact.num.mat <- matrix(0,num.cluster,num.cluster)
	Interact.gene.mat <- matrix(NA,num.cluster,num.cluster)
	Interact.lig.mat <- matrix(NA,num.cluster,num.cluster)
	Interact.rep.mat <- matrix(NA,num.cluster,num.cluster)
	marker.lig.dat <- as.data.frame(matrix(NA,0,ncol(marker.dat)))
	marker.rep.dat <- as.data.frame(matrix(NA,0,ncol(marker.dat)))
	for (lig.idx in 1:num.cluster){
		cluster.l <- cluster.level[lig.idx]
		markers.l <- marker.dat[marker.dat$cluster == cluster.l, ]
		ligands <- markers.l[(markers.l$avg_log2FC > logFC.thre) & (markers.l$p_val_adj < p.thre), 'gene']
		for (rep.idx in 1:num.cluster){
			cluster.r <- cluster.level[rep.idx]
			markers.r <- marker.dat[marker.dat$cluster == cluster.r, ]
			receptors <- markers.r[(markers.r$avg_log2FC>logFC.thre) & (markers.r$p_val_adj<p.thre), 'gene']
			pair.valid <- which((ligs %in% ligands) & (reps %in% receptors))
			if (length(pair.valid) > 0){
				Interact.gene.mat[lig.idx,rep.idx] <- paste(paste(ligs[pair.valid],reps[pair.valid],sep='--'),collapse=';')
				Interact.lig.mat[lig.idx,rep.idx] <- paste(ligs[pair.valid],collapse=';')
				Interact.rep.mat[lig.idx,rep.idx] <- paste(reps[pair.valid],collapse=';')
				marker.lig.dat <- rbind(marker.lig.dat,markers.l[markers.l$gene %in% ligs[pair.valid],])
				marker.rep.dat <- rbind(marker.rep.dat,markers.r[markers.r$gene %in% reps[pair.valid],])
			}
			Interact.num.mat[lig.idx,rep.idx] <- length(pair.valid)
		}
	}
	rownames(Interact.num.mat) <- cluster.level
	colnames(Interact.num.mat) <- cluster.level
	Interact.num.dat <- melt(Interact.num.mat, varnames=c('Cell.From','Cell.To'),value.name="LR.Number",  na.rm=TRUE)

	rownames(Interact.gene.mat) <- cluster.level
	colnames(Interact.gene.mat) <- cluster.level
	Interact.gene.dat <- melt(Interact.gene.mat,varnames=c('Cell.From','Cell.To'),value.name="LR.Info", na.rm = TRUE)
	
	rownames(Interact.lig.mat) <- cluster.level
	colnames(Interact.lig.mat) <- cluster.level
	Interact.lig.dat <- melt(Interact.lig.mat,varnames=c('Cell.From','Cell.To'),value.name="L.Info", na.rm = TRUE)

	rownames(Interact.rep.mat) <- cluster.level
	colnames(Interact.rep.mat) <- cluster.level
	Interact.rep.dat <- melt(Interact.rep.mat,varnames=c('Cell.From','Cell.To'),value.name="R.Info", na.rm = TRUE)
	
	Interact.gene.dat$L.Info <- Interact.lig.dat$L.Info
	Interact.gene.dat$R.Info <- Interact.rep.dat$R.Info
	rownames(Interact.gene.dat) <- 1:nrow(Interact.gene.dat)

	LRpair <- Interact.gene.dat$LR.Info
	names(LRpair) <- paste(Interact.gene.dat$Cell.From,Interact.gene.dat$Cell.To,sep='--')
	lr.split.list <- sapply(LRpair, function(LR){
		strsplit(LR,split=';')
	})
	ident.pair <- names(lr.split.list)
	for (each.pair in ident.pair){
		lr.split.list[[each.pair]] <- paste(each.pair, lr.split.list[[each.pair]],sep='--')
	}
	lr.split.list <- sapply((unlist(lr.split.list)), function(ident.LR.info){
		strsplit(ident.LR.info, split='--')
		})
	lr.unfold.dat <- as.data.frame(t(as.data.frame(lr.split.list)))
	colnames(lr.unfold.dat) <- c('Cell.From','Cell.To','Ligand','Receptor')
	rownames(lr.unfold.dat) <- paste(lr.unfold.dat$Cell.From,lr.unfold.dat$Cell.To,lr.unfold.dat$Ligand,lr.unfold.dat$Receptor, sep='--')

	marker.lig.dat <- unique(marker.lig.dat)
	marker.rep.dat <- unique(marker.rep.dat)
	Interact <- list(InteractNumer=Interact.num.dat,InteractGene=Interact.gene.dat, InteractGeneUnfold=lr.unfold.dat, markerL=marker.lig.dat, markerR=marker.rep.dat, logFC.thre=logFC.thre, p.thre=p.thre, species=species)
	return(Interact)
}

#' To present a dot plot for specific ligand-receptor pairs in specific clusters
#' @param marker.dat Data frame containing information of marker genes
#' @param Interact Interact list returned by findLRpairs
#' @param ligand.ident Vector containing the ligand ident
#' @param receptor.ident Vector containing the receptor ident
#' @return Dotplot showing the ligand-receptor interaction between the selected ligand.ident and receptor.ident
#' @export
dotPlot <- function(marker.dat, Interact, ligand.ident=NULL, receptor.ident=NULL){
	options(stringsAsFactors=F)
	if (is.null(ligand.ident) & is.null(receptor.ident)){
		stop("either ligand.ident or ligand.ident need to be asigned")
	}
	if (length(ligand.ident)>1 & length(receptor.ident)>1){
		stop("specify one cluster for ligand or receptor analysis")
	}

	# get the InteractGene dataframe
	inter.gene.dat <- Interact$InteractGene

	if (length(ligand.ident)==1){
		inter.gene.dat$Xaxis <- inter.gene.dat$Cell.To
		shape <- 16
	}else{
		inter.gene.dat$Xaxis <- inter.gene.dat$Cell.From
		shape <- 17
	}

	logFC.thre <- Interact$logFC.thre
	p.thre <- Interact$p.thre

	inter.ident.dat <- inter.gene.dat
	if (!is.null(receptor.ident)){
		inter.ident.dat <- inter.ident.dat[inter.ident.dat$Cell.To %in% receptor.ident,]
	}
	if (!is.null(ligand.ident)){
		inter.ident.dat <- inter.ident.dat[inter.ident.dat$Cell.From %in% ligand.ident,]
	}

	lr.ident.pair <- inter.ident.dat$LR.Info
	lr.ident.pair <- unlist(sapply(lr.ident.pair,function(x){strsplit(x,split=';')}))
	names(lr.ident.pair) <- NULL
	lr.ident.pair <- unique(lr.ident.pair)

	lr.ident.split.pair <- sapply(lr.ident.pair, function(x){strsplit(x,split='--')})
	ident.ligs <- unlist(lapply(lr.ident.split.pair,function(x){x[1]}))
	ident.reps <- unlist(lapply(lr.ident.split.pair,function(x){x[2]}))

	inter.ident.unfold.dat <- bind_rows(replicate(length(lr.ident.split.pair), inter.ident.dat, simplify = FALSE))
	inter.ident.unfold.dat$Lig <- rep(ident.ligs, each=nrow(inter.ident.dat))
	inter.ident.unfold.dat$Rep <- rep(ident.reps, each=nrow(inter.ident.dat))
	inter.ident.unfold.dat$LR.Info <- paste0(inter.ident.unfold.dat$Lig,' --> ',inter.ident.unfold.dat$Rep)
	### fc.lr and p.lr are to save the measured FC and pval of each LR pair 
	fc.lr <- c()
	p.lr <- c()
	for (each.row in 1:nrow(inter.ident.unfold.dat)){
		current.from <- inter.ident.unfold.dat[each.row,'Cell.From']
		current.to <- inter.ident.unfold.dat[each.row,'Cell.To']
		current.lig <- inter.ident.unfold.dat[each.row,'Lig']
		current.rep <- inter.ident.unfold.dat[each.row,'Rep']

		fc.lig <- marker.dat[marker.dat$cluster==current.from,][current.lig, 'avg_log2FC']
		p.lig <- marker.dat[marker.dat$cluster==current.from,][current.lig, 'p_val_adj']
		
		fc.rep <- marker.dat[marker.dat$cluster==current.to,][current.rep, 'avg_log2FC']
		p.rep <- marker.dat[marker.dat$cluster==current.to,][current.rep,'p_val_adj']

		### if the ligand have a FC > logFC.thre and a p.adj < p.thre
		### and if the receptor have a FC > logFC.thre and a p.adj < p.thre
		if ((!is.na(fc.lig)) & (fc.lig > logFC.thre) & (!is.na(fc.rep)) & (fc.rep > logFC.thre) & (p.lig < p.thre) & (p.rep < p.thre)){
			fc.lr <- c(fc.lr, fc.lig*fc.rep)
			p.lr <- c(p.lr, 1-(1-p.lig)*(1-p.rep))
		}else{
			fc.lr <- c(fc.lr, NA)
			p.lr <- c(p.lr, NA)
		}
		
	}

	inter.ident.unfold.dat$Log2FC_LR <- fc.lr
	inter.ident.unfold.dat$P_LR <- p.lr
	inter.ident.unfold.dat$Log_P_adj <- -log10(p.adjust(p.lr, method='BH'))
	inter.ident.unfold.dat[which(inter.ident.unfold.dat$Log_P_adj > 30), 'Log_P_adj'] <- 30 


	if (length(ligand.ident)==1){
		x.title <- paste0('Receptor clusters for cluster ',ligand.ident)
	}else{
		x.title <- paste0('Ligand clusters to cluster ',receptor.ident)
	}

	plot <- ggplot(inter.ident.unfold.dat,aes(as.character(Xaxis),LR.Info)) + 
		geom_point(aes(size=Log2FC_LR,col=Log_P_adj),shape=shape) +
		scale_colour_gradient(low="green",high="red") + 
		labs(color='-log10(p.adjust)',size='Log2FC',x=x.title,y="")
	return(plot)
}


#' To find those pathways in which the genesets show overlap with the marker ligand and receptor genes in our dataset
#' @param Interact Interact list returned by findLRpairs
#' @param category Character to indicate which pathway to investigate; one of "go", "kegg", 'wikipathway', and "reactome", or "all" for all pathways
#' @return Interact list containing the ligand-receptor interaction information and the pathways showing overlap with the marker ligand and receptor genes in the dataset
#' @export
findLRpath <- function(Interact, category='all'){
	options(stringsAsFactors=F)
	data('cellTalkData', package=('celltalk'))
	if (category=='all'){
		path.list <- cellTalkData$DataPathway[[Interact$species]]
	}else if(category=='wikipathway'){
		path.list <- cellTalkData$DataWikiPathway[[Interact$species]]
	}
	marker.lig.dat <- Interact$markerL
	marker.rep.dat <- Interact$markerR
	lr.gene <- unique(c(marker.lig.dat$gene,marker.rep.dat$gene))
	
	which.overlap.list <- unlist(
		lapply(path.list, function(x){ 
			if(any(x %in% lr.gene)){ return(TRUE) }else{return(FALSE)}
		}))
	path.lr.list <- path.list[which(which.overlap.list)]
	Interact$pathwayLR <- path.lr.list
	return(Interact)
}

#' To find different enriched pathway between two group cells
#' @param Interact Interact list returned by findLRpath
#' @param gsva.mat Matrix containing the pathway enrichment sorces, with rows representing pathways and columns representing cells. Pathway scores are usually computed from gsva, or other methods aiming to measure the pathway enrichment in cells
#' @param ident.lable Vector indicating the identity lables of cells, and the order of lables are required to match order of cells (columns) in the gsva.mat
#' @param select.ident.1 Identity class to define cells for group 1
#' @param select.ident.2 Identity class to define cells for group 2 for comparison; if 'NULL', use all other cells for comparison
#' @param method Method used for differential enrichment analysis, either 't.test' of 'wilcox.test'
#' @return Dataframe including the statistic result comparing the pathway enrichment sorces between group 1 and group 2, the significant recetor and ligand of group 1 in the pathways, and the corresponding up stream identity class which interact with group 1 by releasing specific ligand
#' @export
diffLRpath <- function(Interact,gsva.mat,ident.lable,select.ident.1,select.ident.2=NULL, method='t.test'){
	options(stringsAsFactors=F)

	path.lr.list <- Interact$pathwayLR
	if (is.null(path.lr.list)){
		stop("no pathway detected, run findLRpath befor diffLRpath")
	}
	# marker.lig.dat <- Interact$markerL
	# marker.lig.dat <- marker.lig.dat[marker.lig.dat$cluster==select.ident,]
	marker.ident1.rep.dat <- subset(Interact$InteractGeneUnfold, Cell.To %in% select.ident.1)

	if (nrow(marker.ident1.rep.dat)==0){
		stop("there is no significant receptor in the selected ident")
	}

	### find those pathways in which genesets hava overlap with the marker receptors of the selected cluster
	ident.rep.gene <- marker.ident1.rep.dat$Receptor
	which.overlap.list <-lapply(path.lr.list, function(x){
		overlap.idx <-  x %in% ident.rep.gene
			if(any(overlap.idx)){ 
				return(paste(x[overlap.idx],collapse=',')) 
			}else{
				return(FALSE)
		}
	})
	overlap.rep.list <- which.overlap.list[which(which.overlap.list!='FALSE')]
	### since we set min.sz in the gsva function, there are some pathways not calculated in the gsva process, we shall remove those pathways
	overlap.rep.list <- overlap.rep.list[names(overlap.rep.list) %in% rownames(gsva.mat)]

	### t.test or wilcox.text for the pathways for the selected cluster
	gsva.ident.mat <- gsva.mat[names(overlap.rep.list),]
	group <- as.character(ident.lable)
	test.res.dat <- pathTest(gsva.ident.mat, group, select.ident.1, select.ident.2, method)

	
	### to find the upstream ident and ligand the ident.1 recieved 
	test.res.dat$cell.up <- NA
	test.res.dat$ligand.up <- NA
	test.res.dat$receptor.in.path <- unlist(overlap.rep.list)
	for (each.row in 1:nrow(test.res.dat)){
		each.rep <- test.res.dat[each.row,'receptor.in.path']
		rep.vec <- strsplit(each.rep, split=',')[[1]]
		for (each.rep in rep.vec){
			each.unfold.rep <- subset(marker.ident1.rep.dat, Receptor==each.rep)
			test.res.dat[each.row,'cell.up'] <- paste(each.unfold.rep$Cell.From, collapse=';')
			if (length(unique(each.unfold.rep$Ligand)) > 1){
				test.res.dat[each.row,'ligand.up'] <- paste(each.unfold.rep$Ligand, collapse=';')
			}else{
				test.res.dat[each.row,'ligand.up'] <- each.unfold.rep$Ligand[1]
			}
		}
	}

	### to find the ligand in the same pathway
	marker.lig.dat <- Interact$markerL
	ident.lig.vec <- marker.lig.dat[marker.lig.dat$cluster %in% select.ident.1,'gene']
	### for each pathway in the DEG result, find which pathway show overlap with marker ligands of the selected ident
	ident.lig.in.path <- sapply(test.res.dat$description, function(eachPath){
		each.set <- path.lr.list[[eachPath]]
		if (any(ident.lig.vec %in% each.set)){
			return(paste(ident.lig.vec[ident.lig.vec %in% each.set],collapse=','))
		}else{
			return(NA)
		}
	})
	test.res.dat$ligand.in.path <- unlist(ident.lig.in.path)

	return(test.res.dat)
}

#' To find the downstream identity class of specific ligand released by specific upstream identity class
#' @param Interact Interact list returned by findLRpairs
#' @param select.ident Upstream identity class; if 'NULL', use all identity classes
#' @param select.ligand Ligand released by upstream identity class; if 'NULL', use all ligands that are markers for the selected upstream identity class
#' @return Dataframe including the interaction information
#' @export
findReceptor <- function(Interact, select.ident=NULL, select.ligand=NULL){
	options(stringsAsFactors=F)

	if (is.null(select.ident) & is.null(select.ligand)){
		stop("either a select.ident or a select.ligand need to be asigned")
	}
	if (is.null(select.ligand)){
		ident.down.dat <- subset(Interact$InteractGeneUnfold, Cell.From==select.ident)
	}else if(is.null(select.ident)){
		ident.down.dat <- subset(Interact$InteractGeneUnfold, Ligand==select.ligand)
	}else{
		ident.down.dat <- subset(Interact$InteractGeneUnfold, Cell.From==select.ident & Ligand==select.ligand)
	}

	if (nrow(ident.down.dat)==0){
		stop("no downstream ident found for the selected ident and ligand")
	}
	return(ident.down.dat)
}

#' Differential enrichment analysis by t.test or wilcox.txt
#' @param gsva.ident.mat Matrix of pathway scores, pathway * cell
#' @param group Vector of group labels of cells
#' @param select.ident.1 Identity class 1
#' @param select.ident.2 Identity class 2 for comparison
#' @param method Method for hypothesis test, either 't.test' or 'wilcox.test'
#' @return Dataframe including the statistic result
#' @export
pathTest <- function(gsva.ident.mat, group, select.ident.1, select.ident.2=NULL, method='t.test'){
	if (method=='t.test'){
		if (is.null(select.ident.2)){
		t.result <- apply(gsva.ident.mat,1,function(geneExpr){
			t.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[!(group %in% select.ident.1)])}
		)
		}else{
			t.result <- apply(gsva.ident.mat,1,function(geneExpr){
				t.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[group %in% select.ident.2])}
			)
		}
		# t = testRes$statistic, df = testRes$parameter, mean.1 = testRes$estimate[1], mean.2 = testRes$estimate[2], pVal = testRes$p.value
		test.res.dat <- as.data.frame(lapply(t.result,function(testRes){
			c(testRes$estimate[1]-testRes$estimate[2],testRes$estimate[1],testRes$estimate[2],testRes$statistic,testRes$parameter,testRes$p.value)
		}))
		test.res.dat <- as.data.frame(t(test.res.dat))
		colnames(test.res.dat) <- c('mean.diff','mean.1','mean.2','t','df','p.val')
	}else if(method=='wilcox.test'){
		if (is.null(select.ident.2)){
			wil.result <- apply(gsva.ident.mat,1,function(geneExpr){
			wilcox.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[!(group %in% select.ident.1)])})

			wil.median <- apply(gsva.ident.mat, 1, function(geneExpr){
				median.1 <- median(geneExpr[group %in% select.ident.1])
				median.2 <- median(geneExpr[!(group %in% select.ident.1)])
				median.diff <- median.1 - median.2
				c(median.diff, median.1, median.2)
			})
			wil.median <- as.data.frame(t(wil.median))
			colnames(wil.median) <- c('median.diff','median.1','median.2')

		}else{
			wil.result <- apply(gsva.ident.mat,1,function(geneExpr){
				wilcox.test(x=geneExpr[group %in% select.ident.1],y=geneExpr[group %in% select.ident.2])})
			wil.median <- apply(gsva.ident.mat, 1, function(geneExpr){
				median.1 <- median(geneExpr[group %in% select.ident.1])
				median.2 <- median(geneExpr[group %in% select.ident.2])
				median.diff <- median.1 - median.2
				c(median.diff, median.1, median.2)
			})
			wil.median <- as.data.frame(t(wil.median))
			colnames(wil.median) <- c('median.diff','median.1','median.2')
		}
		# W = testRes$statistic, pVal = testRes$p.value
		test.res.dat <- as.data.frame(lapply(wil.result,function(testRes){
			c(testRes$statistic,testRes$p.value)
		}))
		test.res.dat <- as.data.frame(t(test.res.dat))
		colnames(test.res.dat) <- c('W','p.val')

		test.res.dat <- cbind(wil.median,test.res.dat)
	}else{
		stop("select t.test or wilcox.test to conduct differential analysis")
	}
	
	test.res.dat$p.val.adj <- p.adjust(test.res.dat$p.val, method='BH')
	test.res.dat$description <- rownames(test.res.dat)
	return(test.res.dat)
}