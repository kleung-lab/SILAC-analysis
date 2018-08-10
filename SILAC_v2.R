# script to calculate log2 ratios and p values from skyline export
# July 7,2017 Kevin Leung
# August 9, 2018 updated
# merged Cal_ratio, volplot, and box_temp scripts into one script
# Functions: calRatio3, log2dist, and volplot, HLpep, HLpro
# Dependancies - ggplot2, reshape2, ggrepel
# 

### calRatio3
# reads in skyline report files, and analyze data according to a parameter file
# skyline reports can either be in a folder, as a single file, or an imported R object
# parameter file defines each experiment in the skyline report (searched by raw file name)

calRatio3 <- function (input = " ", para = " ", analysis_name = " ", extra = F){
  # set up empty dataframe objects
  skyraw <- data.frame()
  skypep <- data.frame()
  skypro <- data.frame()
  skypro_ind <- data.frame()
  
  # read input files
  # if input ends in /, it will read and concatenate the raw files to skyraw. 
  # if input ends in .csv, it will read the csv file to skyraw
  # if input is a data object, it will copy it to skyraw
  if (grepl("/$",input)==T){ 
    f <- list.files(input, pattern = ".csv")
    for (i in f){
      skyraw <- rbind(skyraw, read.csv(file = paste(input,i, sep=""), stringsAsFactors=FALSE))
    }
    cat("Input is a directory, the following skyline files are loaded -", f, sep="\n")
  } else if (grepl(".csv$",input)==T){
    skyraw <- read.csv(file = input, stringsAsFactors=FALSE)
  } else skyraw <- input 
  
  # #read in extracellularome file
  # # omitted for general purpose uses
  # extra_lib <- read.csv("~/Box Sync/MS/MSpeaklist/lib/human_extra.csv", stringsAsFactors = F)
  # extra_lib <- rbind(extra_lib,read.csv("~/Box Sync/MS/MSpeaklist/lib/mus_extra.csv", stringsAsFactors = F) )
  # if(extra == T){
  # skyraw <- skyraw[skyraw$Protein.Accession %in% extra_lib$Acc,]
  # }
  
  
  # read the number of unique experiments from skyline report
  m <- as.character(unique(skyraw[["File.Name"]]))
  if (length(m) == 0){cat("Input is invalid")}
  
  # read parameter file
  #asks for a version number
  pm <- read.csv(file = para, stringsAsFactors=FALSE)
  # if (is.null(analysis_name) == T){
  # analysis_name <- readline (prompt = "analysis name (unique): ")
  # }
  # calculate log2 l/h ratio - filter things out later
  skyraw$lh <- as.numeric(skyraw$light.Total.Area)/as.numeric(skyraw$heavy.Total.Area)
  skyraw[is.nan(skyraw$lh), "lh"] <- NA
  skyraw[is.infinite(skyraw$lh), "lh"] <- NA # set inf to NA
  skyraw[!is.na(skyraw$lh) & skyraw$lh == 0, "lh"] <- NA # set 0 to NA
  skyraw$lh <- log2(skyraw$lh) #transform into log2 ratio
  
  # for each dataset found in skyline report, specificed by m -
  # only process data if it's going to be used (parameter file)
  # flip light/heavy ratio if the experiment is H (parameter file)
  # do l/h filtering according to dotp value (datafile)
  # filter out deamidated peptides if it's PNG (parameter file)
  
  for (i in 1:length(m)){
    # (parameter "used") search for unique replicate names in parameter raw
    if (pm[grepl(m[i],pm[["raw"]]),"used"] == T) {
      # extract parameter for this file
      s.pm <- pm[grepl(m[i],pm[["raw"]]),]
      # extract dataset for this file
      s.data <- skyraw[grepl(m[i], skyraw[["File.Name"]]),]
      # append with unique experiment name and grouped exp_ref
      s.data$exp_name <- s.pm$exp_name
      s.data$exp_ref <- s.pm$exp_ref
      
      ### FILTERING based on dotp, PNGase fraction, and minimum number of peptides
      # data for a given dotp, can filter based on either H or L ratio > dotp or BOTH > dotp 
      # make h and l filters as logical table based on dotp
      l <- s.data[["light.Isotope.Dot.Product"]] > s.pm$dotp
      h <- s.data[["heavy.Isotope.Dot.Product"]] > s.pm$dotp
      # make fil logical filter based on dotp.fil
      # if dotp.fil is OR, either heavy or light dotp has to be greater than dotp 
      # if dotp.fil is blank, BOTH heavy and light dotp has to be greater than specified
      if (s.pm$dotp.fil == "OR"){
          fil.dotp <- l|h
      } else if (s.pm$dotp.fil == "AND"){
          fil.dotp <- l&h
      }

      # filter out non-deamidated peptides for PNG fractions
      # using parameter "fraction" to make fil.PNG table
      if (s.pm$fraction == "PNG"){
        fil.PNG <- grepl("N.+1.", s.data[["Peptide.Modified.Sequence"]])
      }else {
          fil.PNG <- rep(T, nrow(s.data))
      }
      # filtering out proteins based on number of well quantified peptides 
      # use parameter "peptides", can be set separately for PNGase and tryptic fraction
      # names of proteins where enrich.fil > specified
      fil.pro <- s.data[fil.dotp & fil.PNG,"Protein.Accession"] #well quantified proteins based on fil.dotp and fil.PNG
      fil.pro <- table(fil.pro) >= s.pm$peptides #table of proteins with well quantified peptides greater than specified
      fil.pro <- s.data$Protein.Accession %in% names(fil.pro)[fil.pro] #subset of proteins that matches those in fil.pro
      
      # merging the three filters together as fil and append it to s.data
      fil <- fil.pro & fil.dotp & fil.PNG
      s.data$fil <- fil
      
      ###### centering and calculating enrichment ratio
      # centering data (sd and center)
      # calculate mean and sd based on filtered peptides
      lh.mean <- mean(s.data[fil,"lh"], na.rm=T)
      lh.sd <- sd(s.data[fil,"lh"], na.rm=T)
      # calculate lh.norm based on center and sd parameter
      if (s.pm$center == T && s.pm$sd ==T){
          s.data$lh.norm <- (s.data$lh-lh.mean)/lh.sd
      } else if (s.pm$center == T && s.pm$sd ==F){
          s.data$lh.norm <- s.data$lh-lh.mean
      } else if (s.pm$center == F && s.pm$sd ==F){
          s.data$lh.norm <- s.data$lh
      }
      # (parameter "exp_label")
      # flip l/h to h/l if exp_label is H
      if (s.pm$exp_label == "H"){
          s.data$enrich <- -s.data$lh.norm
      } else if (s.pm$exp_label == "L"){
       s.data$enrich <- s.data$lh.norm}

      ######
      # Calculate median ratio and p value using only data that're filtered
      # p value test is mann whitney one sample test between all the values for a single protein compared to 0 - this is the correct test for median ratio statistics test. t-test looks at mean value between two groups of values and is not suitable for SILAC ratio - worthwhile thinking through this. 
      # aggregate "enrich" based on "Protein.Accession" and "Protein.Gene", 
      # data is s.data[fil,] 
      # calculates median, p value, and replicate count
      rep_ratio <- as.data.frame(as.list(suppressWarnings(aggregate(enrich ~ Protein.Accession + Protein.Gene, s.data[fil,], function(x) c(rep.med = median (x), rep.p = -log10(wilcox.test(x, mu=0)$p.value), rep.count = length(x))))))
      # writing the protein level data as skypro_ind for each replicate, useful for comparing individual replicate
      s.pro <- rep_ratio
      s.pro$exp_name <- s.pm$exp_name
      s.pro$exp_ref <- s.pm$exp_ref
      s.pro$raw_id <- s.pm$raw
      skypro_ind <- rbind(skypro_ind,s.pro)
      
      # appends replicate protein data to the peptide level "s.data" by protein accession number
      rep_ratio$Protein.Gene <- NULL
      s.data <- merge(s.data, rep_ratio, by="Protein.Accession", all.x = T)
      
      #concatenate all data into skypep from s.data
      skypep <- rbind(skypep,s.data) 

    }
  }
  
  # Calculate median ratio and p value at the experiment level
  # same function as the one at replicate level, except now it aggregates by exp_ref in addition to protein accession and protein gene name
  # calculates median, p value, and replicate count
    skypro <- as.data.frame(as.list(suppressWarnings(aggregate(enrich ~ Protein.Accession +Protein.Gene+ exp_ref, skypep[skypep$fil,], function(x) c(exp.med = median (x), exp.p = -log10(wilcox.test(x, mu=0)$p.value), exp.count = length(x))))))
   
  # rename columns of file for exporting into sql database
  # May be easier to do this in separate function??
  # sky2sql <- read.csv("./sky2sql.csv", stringsAsFactors = F)
  #  
  # for (i in 1:nrow(sky2sql)){
  #   names(skypep)[names(skypep) == sky2sql[i,"sky"]] <- sky2sql[i,"sql"]
  #   names(skypro)[names(skypro) == sky2sql[i,"sky"]] <- sky2sql[i,"sql"]
  #   names(skypro_ind)[names(skypro_ind) == sky2sql[i,"sky"]] <- sky2sql[i,"sql"]
  #   # class(skypep)[names(skypep) == sky2sql[i,"sql"]] <- sky2sql[i,"class"]
  # }
  
    
    #### append each file with given analysis name, not sure if it's useful, maybe to keep track of different parameters
  pm[,"analysis"] <- analysis_name
  skypep[,"analysis"] <- analysis_name
  skypep <- skypep[order(skypep$exp_name,skypep$Protein.Gene),]
  skypro[,"analysis"] <- analysis_name
  skypro_ind[,"analysis"] <- analysis_name

  # write files to global environment
  assign("parameter", pm, .GlobalEnv)
  assign("skypep", skypep, .GlobalEnv)
  assign("skypro",skypro,.GlobalEnv)
  assign("skypro_ind",skypro_ind,.GlobalEnv)
}

# Function to plot log2 ratio distribution (before normalization, using lh from skypep)
# box style or histograme style
log2dist <- function (skypep, style = "box"){
    if (style == "box"){
        par(mfrow=c(1,1))
        d <- skypep[skypep[["fil"]]==T,]
        d <- d[d$lh<10,]
        d <- d[d$lh>-10,]
        d <- d[!is.na(d$lh),]
        par(mfrow=c(1,1))
        
        library(ggplot2)
        # library(RColorBrewer)
        
        p <- ggplot(d,aes(x=exp_name, y=lh, fill = exp_ref))+ geom_boxplot()+theme_classic()
        p + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
            geom_hline(yintercept = c(-1,1),colour="grey", linetype = "longdash") 
        
    
        
    } else if (style == "hist"){
        par(mfrow=c(2,2))
        
        # plotting log2 ratios distribution (before normalization, only showing the filtered peptides)
        for (i in unique(skypep$exp_name)){
            d <- skypep[skypep[["exp_name"]]==i,]
            d <- d[d[["fil"]]==T,]
            lh.mean <- mean(d[,"lh"], na.rm=T)
            lh.sd <- sd(d[,"lh"], na.rm=T)
            hist(d[,"lh"], xlab = "log2 l/h ratio",ylab="peptide number", main=i, xlim=c(-8,8), breaks = c(-20,seq(-10,10,1),20), freq=T, xaxt="n")
            axis(side=1, at=seq(-8,8,by=2))
            abline(v=lh.mean, col="red")
            legend("topright", paste("mean =",round(lh.mean,2), "\n  SD =",round(lh.sd,2)), bty="n")  
                
        }

    }
}
# Plotting light vs heavy enrichemnt ratios
HLpro <- function (skypep) {
    library(reshape2)
    for (i in unique(skypep$exp_ref)){
        d <- skypep[skypep[["exp_ref"]]==i,]
        d <- d[d[["fil"]]==T,]
        pro_comp <- dcast(d, Protein.Gene ~ exp_name, value.var = "enrich.rep.med", mean)
        pro_comp$Protein.Gene <- NULL
        pro_comp[pro_comp > 5 | pro_comp< -5] <- NA
        
        pairs(pro_comp, upper.panel = function(x,y){
            points(x,y)
            abline(lm(y~x), col='red')
        },
        lower.panel = function(x,y) {
            legend("top", xjust = 0.5, bty = "n", legend = paste("fit R2 =\n",format(summary(lm(y~x))$adj.r.squared, digits=3), "\nPearson r =\n",format(cor(x,y, use = "complete.obs", method = "pearson"), digits = 3)), cex =1.2, pt.cex = 1)
            
            
        }, main = paste("Replicate comparison of protein enrichment level in", i)
        )
        
    }
}
HLpep <- function (skypep) {
    library(reshape2)
    for (i in unique(skypep$exp_ref)){
        d <- skypep[skypep[["exp_ref"]]==i,]
        d <- d[d[["fil"]]==T,]
        pep_comp <- dcast(d, Peptide.Modified.Sequence ~ exp_name, value.var = "enrich", mean)
        pep_comp$Peptide.Modified.Sequence <- NULL
        pep_comp[pep_comp > 5 | pep_comp< -5] <- NA
        
        pairs(pep_comp, upper.panel = function(x,y){
            points(x,y)
            abline(lm(y~x), col='red')
        },
        lower.panel = function(x,y) {
            legend("center", bty = "n", legend = paste("R2 =",format(summary(lm(y~x))$adj.r.squared, digits=4)))
            
            
        }, main = paste("Replicate comparison of peptide enrichment level in", i)
        )
        
    }
}



# plot volcano plot by experiment (based on exp_rep)
Volplot <- function(input, custom = " ", r.cutoff = 1, p.cutoff=0.05, customlab = "normal"){
    f <- input #data with three columns – gene, ratio, and p
    colnames(f) <- c("acc", "gene", "exp_ref", "ratio", "p", "analysis")
    f$exp_ref <- as.character(f$exp_ref)
    
    
    library(ggplot2)
    library(ggrepel)
    set.seed(42)
    l <- unique(f$exp_ref)
    for (i in 1:length(l)){
        g <- NULL
        g <- f[grepl(l[i], f$exp_ref),]
        
        #adding a column for coloring and labeling scheme using nested ifelse statments
        g$Significant <- ifelse(as.character(g$gene) %in% custom & g$p > -log10(p.cutoff) & (g$ratio < -r.cutoff| g$ratio > r.cutoff),"4",
                                ifelse (as.character(g$gene) %in% custom,"3", 
                                         ifelse(g$p > -log10(p.cutoff) & g$ratio < -r.cutoff, "2", 
                                                ifelse(g$p > -log10(p.cutoff) & g$ratio>r.cutoff, "1", "0"))) )
        
        # pdf(paste("./_test/vol_",cells[i],"_v2.pdf",sep=""), 8,8, onefile=F)
        p <- ggplot(g,aes(ratio,p))+
            geom_point(aes(color = Significant)) +
            scale_color_manual(values = c("grey", "green","red","blue","purple")) +
            geom_hline(yintercept = -log10(p.cutoff),colour="grey", linetype = "longdash")+
            geom_vline(xintercept = c(r.cutoff,-r.cutoff),colour="grey", linetype = "longdash")+
            theme_bw(base_size = 18) +
            theme(legend.position = "none")+
            labs(x = expression(paste("median log2 enrichment ratio")), y = "-log10 p values")+
            ggtitle(paste("SILAC ratio of", l[i]))+
            
            theme(plot.title = element_text(hjust = 0.5),
                  panel.border = element_blank(),
                  panel.background = element_blank(),
                  panel.grid.major.x = element_blank(),
                  panel.grid.minor.x = element_blank(),
                  panel.grid.major.y = element_blank(),
                  panel.grid.minor.y = element_blank(),
                  axis.line.x = element_line(),
                  axis.line.y= element_line())+
            annotate("text", label = paste("p value <" , p.cutoff,
                                           "\nratio cutoff =",r.cutoff,
                                           # "\ndotp >", dotp,
                                           "\ntotal proteins =",nrow(g),
                                           "\nratio SD =", round(sd(g$ratio),2)),
                     # "\ntotal rep =",ncol(z.vol.temp),
                     # "\naverage peptides per rep =",round(nrow(z.vol[!is.na(z.vol$ratio),])/ncol(z.vol.temp))),
                     x = max(g[g$p > -log10(p.cutoff),"ratio"])*0.8, y = max(g$p)*0.9) 
        # cut off data points with p value below threshold
        p <- p + xlim(min(g[g$p > -log10(p.cutoff),"ratio"]-1),max(g[g$p > -log10(p.cutoff),"ratio"])+1)
        
        if(customlab == "normal"){ 
            q <- p + geom_text_repel(
                data = subset(g, Significant ==1),
                aes(label = gene),
                nudge_x = 0.1,
                nudge_y = 0.1
            ) +
                geom_text_repel(
                    data = subset(g, Significant ==2),
                    aes(label = gene),
                    nudge_x = -0.1,
                    nudge_y = 0.1 )
        } else if (customlab == "custom"){ 
            q <- p+ geom_text_repel(
                data = g[g$Significant >= 3, ],
                aes(label = gene), fontface="bold") 
        } else if (customlab == "orfeome"){
            q <- p+ geom_text_repel(
                data = subset(g, Significant ==4),
                aes(label = gene), fontface="bold") 
        }
        print(q)
        assign(paste(deparse(substitute(input)),"-",l[i], "-graph",sep=""),q, envir=globalenv())
    }
}





