# How to run the script:
#Rscript alpha_diversity.R diversity_check.csv "Treatment 1" "Treatment 2" "Treatment 3" "Treatment 4" treatment "Clinical Outcome" testing_kw.pdf kw

start.time <- Sys.time()
require('ggplot2')

require('grid')
require('gridExtra')


# Multiple plot function: (http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/)
multiplot <- function(..., plotlist=NULL, file, cols=2, layout=NULL) {
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  numPlots = length(plots)
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                    ncol = cols, nrow = ceiling(numPlots/cols))
  }
 if (numPlots==1) {
    print(plots[[1]])
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}



# Function to re-order a dataframe based on a user-defined order of a specific column entries
dataframe_filter <- function(namearrange,meta,divint,MYdata){
	rearrlist <- lapply(namearrange,function(x) eval(parse(text = paste("subset(MYdata,", meta, "==x,select=c(",meta,",",divint,"))"))))
	newdataf <- do.call("rbind", rearrlist)
	return (newdataf)	
}


# Function to create dotplots with ANOVA and Kruskal-Wallis p-values
dot_plotting <- function(MYdata,namearrange,meta,ylabel,xlabel,comp1,comp2,comp3,divint,maxval,sigtest){
	MYdata[,meta] = factor(MYdata[,meta],namearrange)
	comp1df <- dataframe_filter(comp1,meta,divint,MYdata)
	comp2df <- dataframe_filter(comp2,meta,divint,MYdata)
	comp3df <- dataframe_filter(comp3,meta,divint,MYdata)
	
	pval.kw1 <- signif(kruskal.test(as.formula(paste(divint,"~",meta)), data=comp1df)[3][[1]],3)
	pval.kw2 <- signif(kruskal.test(as.formula(paste(divint,"~",meta)), data=comp2df)[3][[1]],3)
	pval.kw3 <- signif(kruskal.test(as.formula(paste(divint,"~",meta)), data=comp3df)[3][[1]],3)
	pval.anov1 <- signif(summary(aov(as.formula(paste(divint,"~",meta)), data=comp1df))[[1]][[5]][1],3)
	pval.anov2 <- signif(summary(aov(as.formula(paste(divint,"~",meta)), data=comp2df))[[1]][[5]][1],3)
	pval.anov3 <- signif(summary(aov(as.formula(paste(divint,"~",meta)), data=comp3df))[[1]][[5]][1],3)
	pval.kw.test1 <- paste("p-value = ",ifelse(pval.kw1 < 0.05,paste(pval.kw1,"*",sep=""),pval.kw1),sep="")		# KW p-value group 1 vs group 4
	pval.kw.test2 <- paste("p-value = ",ifelse(pval.kw2 < 0.05,paste(pval.kw2,"*",sep=""),pval.kw2),sep="")		# KW p-value group 2 vs group 4
	pval.kw.test3 <- paste("p-value = ",ifelse(pval.kw3 < 0.05,paste(pval.kw3,"*",sep=""),pval.kw3),sep="")		# KW p-value group 3 vs group 4
	pval.anova.test1 <- paste("p-value = ",ifelse(pval.anov1 < 0.05,paste(pval.anov1,"*",sep=""),pval.anov1),sep="")		# Anova p-value group 1 vs group 4
	pval.anova.test2 <- paste("p-value = ",ifelse(pval.anov2 < 0.05,paste(pval.anov2,"*",sep=""),pval.anov2),sep="")		# Anova p-value group 2 vs group 4
	pval.anova.test3 <- paste("p-value = ",ifelse(pval.anov3 < 0.05,paste(pval.anov3,"*",sep=""),pval.anov3),sep="")		# Anova p-value group 3 vs group 4
	
	cmd0 <- paste('ggplot(MYdata, aes(x = factor(',meta,'), y =',divint,', fill=factor(',meta,'))) + geom_dotplot(binaxis = "y", stackdir = "center",position="dodge")',sep="")
	cmd1 <- paste(cmd0,' + geom_segment(aes(x = 1, y = ',maxval/0.99,', xend = 2, yend = ',maxval/0.99,'))',sep="")
	cmd2 <- paste(cmd1,' + geom_segment(aes(x = 1, y = ',maxval/0.95,', xend = 3, yend = ',maxval/0.95,'))',sep="")
	cmd3 <- paste(cmd2,' + geom_segment(aes(x = 1, y = ',maxval/0.91,', xend = 4, yend = ',maxval/0.91,'))',sep="")
	cmd4 <- paste(cmd3,' + xlab("',xlabel,'") + ylab("',ylabel,'") + guides(fill=guide_legend(title="Clinical Outcome"))',sep="")
	
	# Code for adjusting axes titles, legend names, etc.
	#theme(plot.title = element_text(lineheight=.8, face="bold",size=30),axis.text=element_text(size=20),
	#axis.title=element_text(size=40,face="bold"),legend.key.size = unit(1, "cm"),legend.text=element_text(size=30),
	#legend.title=element_text(size=30)) + 
	
	col_kw <- ifelse(c(pval.kw1,pval.kw2,pval.kw3) < 0.05,"red","black")
	col_anov <- ifelse(c(pval.anov1,pval.anov2,pval.anov3) < 0.05,"red","black")
	color_anova <- paste('c("',paste(col_anov, collapse = '","'),'")',sep="")
	color_kw <- paste('c("',paste(col_kw, collapse = '","'),'")',sep="")
	
	cmd_kw <- paste(cmd4,' + annotate("text", x = c(1.5,2,2.5), y = c(',maxval/0.97,',',maxval/0.93,',',maxval/0.89,'), label = c("',pval.kw.test1,'","',pval.kw.test2,'","',pval.kw.test3,'"), colour=',color_kw,', size=3)',sep="")
	cmd_anova <- paste(cmd4,' + annotate("text", x = c(1.5,2,2.5), y = c(',maxval/0.97,',',maxval/0.93,',',maxval/0.89,'), label = c("',pval.anova.test1,'","',pval.anova.test2,'","',pval.anova.test3,'"), colour=',color_anova,', size=3)',sep="")
	if (sigtest == "anova"){
		finalfig <- eval(parse(text=cmd_anova))
		return (finalfig)
	} else {
		finalfig <- eval(parse(text=cmd_kw))
		return (finalfig)
	}	
}


# Function to call all other functions
main_function <- function(infile,namearrange,meta,xlabel,comp1,comp2,comp3,diversity_measures,ylabels,outfile,sigtest){
	MYdata <- read.csv(infile, header = T, sep = ",", check.names = T, row.names = 1)
	maxval1 <- ceiling(max(MYdata[,diversity_measures[1]]))
	fg1 <- dot_plotting(MYdata,namearrange,meta,ylabels[1],xlabel,comp1,comp2,comp3,diversity_measures[1],maxval1,sigtest)
	maxval2 <- ceiling(max(MYdata[,diversity_measures[2]]))
	fg2 <- dot_plotting(MYdata,namearrange,meta,ylabels[2],xlabel,comp1,comp2,comp3,diversity_measures[2],maxval2,sigtest)
	maxval3 <- max(MYdata[,diversity_measures[3]])
	fg3 <- dot_plotting(MYdata,namearrange,meta,ylabels[3],xlabel,comp1,comp2,comp3,diversity_measures[3],maxval3,sigtest)
	maxval4 <- ceiling(max(MYdata[,diversity_measures[4]]))
	fg4 <- dot_plotting(MYdata,namearrange,meta,ylabels[4],xlabel,comp1,comp2,comp3,diversity_measures[4],maxval4,sigtest)
	maxval5 <- ceiling(max(MYdata[,diversity_measures[5]]))
	fg5 <- dot_plotting(MYdata,namearrange,meta,ylabels[5],xlabel,comp1,comp2,comp3,diversity_measures[5],maxval5,sigtest)
	pdf(outfile,width=30,height=13)			# width 20,30 , height 8,13
	multiplot(fg1,fg2,fg3,fg4,fg5)
	dev.off()
}


# Command arguments
argv <- commandArgs(TRUE)
infile <- argv[1]
trt1 <- argv[2]
trt2 <- argv[3]
trt3 <- argv[4]
trt4 <- argv[5]
namearrange <- c(trt1,trt2,trt3,trt4)
meta <- argv[6]
xlabel <- argv[7]
comp1 <- c(trt1,trt2)
comp2 <- c(trt1,trt3)
comp3 <- c(trt1,trt4)
outfile <- argv[8]
sigtest <- argv[9]

# Variables not changing
diversity_measures <- c("richness","shannon","evenness","inv_simpson","PD_whole_tree")
ylabels <- c("Richness","Shannon's Diversity","Evenness","Inverse Simpson","Faith's Phylogenetic Diversity")
main_function(infile,namearrange,meta,xlabel,comp1,comp2,comp3,diversity_measures,ylabels,outfile,sigtest)
print (Sys.time() - start.time)
