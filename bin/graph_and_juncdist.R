library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

myparse = function(file,filedesc) {
   dfall = read.table(file,sep=",",header=T)
   
   dfall$sampleID = as.character(dfall$sampleID)
   dfall$sample = as.character(dfall$sample)
   dfall$treat = as.character(dfall$treat)
   dfall$isotype = as.character(dfall$isotype)
   if (length(grep("^type$",colnames(dfall))) > 0) {
      dfall$type = as.character(dfall$type)
      df = dfall
      df = df[grep("mis_[A-Z]..N",df$type,perl=T,invert=T),]
      dm1 = df
   } else {
      dfall$strand = as.character(dfall$strand)
   }
   #    else if (length(grep("^strand2$",colnames(dfall))) > 0) {
   #    dfall$strand2 = as.character(dfall$strand2)
   #    dm1 = dfall
   # } else {
   #    dfall$strand2 = "-"
   #    dfall$strand2 = as.character(dfall$strand2)
   #    dm1 = dfall
   # }
   dm1$myfill = paste(dm1$treat,dm1$isotype)#,"\nN =",dm1$total_read)
   dm1$myfill = as.character(dm1$myfill)
   dm1$myfill = gsub("^([0-9]+)\\. (.+) \\(([AB]\\. .+)\\) [0-9][0-9]\\. (.+)","TREAT is \\1 (\\2)
\\3
BAIT at IgM
PREY at \\4",dm1$myfill,perl=T)
   dm1$file = filedesc
   
   return(dm1)
}


mygraph = function(dm1,dm1desc) {
   mysampletot = length(unique(as.character(dm1$myfill))) * length(unique(as.character(dm1$sample)))
   
   print(paste(dm1desc,", total sample = ",mysampletot,sep=""))
   themebw = theme_bw() + theme(panel.grid=element_blank(),axis.ticks.x=element_blank(),axis.text.x=element_blank())
   themebwwithx = theme_bw() + theme(panel.grid=element_blank())
   
   p1.dm = dm1[grep("^(mis|mh|ins|del)$",dm1$type,perl=T),]
   
   p1.dm$num = gsub("^([0-9]+)\\..+$","\\1",p1.dm$sample,perl=T)
   p1.dm$num = as.numeric(as.character(p1.dm$num))
   
   p1 = ggplot(p1.dm,aes(num,perc_atleast1)) +
      # stat_summary(geom="rect",fun.xmin=num-0.5,fun.xmax=num+0.5,fun.ymin=0,fun.ymax=perc_atleast1,aes(fill=sample),color="black") +
      #geom_rect(aes(xmin=num-0.5,xmax=num+0.5,ymin=0,ymax=perc_atleast1,fill=sample),color="black") +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9),color="black") +
      # stat_summary(geom="bar",fun.y=mean,aes(fill=sample),color="black",stat="identity",position=position_dodge(width=0.9)) +
      facet_grid(type~myfill,scales = "free_y") +
      geom_text(aes(y=perc_atleast1,group=sample,label=perc_atleast1),vjust=0) +
      themebw + xlab("") + ylab("% of read with at least 1 [MIS/INS/DEL/MH]") + ggtitle("% of read with at least 1 [MIS/INS/DEL/MH] of relevant isotypes")#COUNT[i=0->i=n]{mutation >= 1}/COUNT[i=0->i=n]{n} => <Total number of read with at least 1 [MIS_X->Y]>) / <Total read>")
   if (length(grep("combine",dm1desc)) > 0) {
      p1 = p1 + geom_errorbar(aes(group=sample,ymin=perc_atleast1-perc_atleast1.se,ymax=perc_atleast1+perc_atleast1.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=perc_atleast1.pchar))
   }   
   
   p1tot = ggplot(p1.dm[p1.dm$type == "del",],aes(num,total_read)) +
      geom_bar(aes(fill=sample),color="black",stat="identity",position=position_dodge(width=0.9)) +
      facet_grid(.~myfill,scales="free_y") +
      geom_text(aes(y=0.8*total_read,group=sample,label=total_read),stat="identity",position=position_dodge(width=0.9),vjust=0) +
      themebw + xlab("") + ylab("Number of read") + ggtitle("Total read")
   if (length(grep("combine",dm1desc)) > 0) {
      p1tot = p1tot + geom_errorbar(aes(group=sample,ymin=total_read-total_read.se,ymax=total_read+total_read.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=total_read.pchar))
   }   
   
   
   p1.perbp_allread = ggplot(p1.dm,aes(num,mean_perbp_allread)) +
      #geom_rect(aes(xmin=num-0.5,xmax=num+0.5,ymin=0,ymax=mean_perbp_allread,fill=sample),color="black") +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9),color="black") +
      facet_grid(type~myfill,scales = "free_y") +
      geom_text(aes(y=mean_perbp_allread,group=sample,label=mean_perbp_allread),vjust=0) +
      themebw + xlab("") + ylab("mutation/bp") + ggtitle("Stat of [MIS/INS/DEL/MH] of relevant isotypes; mean(all reads)/mean(all bp)\nSUM[i=0->i=n]{mutation}/SUM[i=0->i=n]{bp} => <Total mutation of all reads> / <Total bp of all reads> => ( (1+4+6+2) / N ) / ( (70+65+60+90) / N )")

   if (length(grep("combine",dm1desc)) > 0) {
      p1.perbp_allread = p1.perbp_allread + geom_errorbar(aes(group=sample,ymin=mean_perbp_allread-mean_perbp_allread.se,ymax=mean_perbp_allread+mean_perbp_allread.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=mean_perbp_allread.pchar))
   }   

      p1.perbp_eachread = ggplot(p1.dm,aes(num,1000*mean_perbp_eachread)) +
      #geom_rect(aes(xmin=num-0.5,xmax=num+0.5,ymin=0,ymax=1000*mean_perbp_eachread,fill=sample),color="black") +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9),color="black") +
      facet_grid(type~myfill,scales = "free_y") +
      geom_text(aes(y=1000*mean_perbp_eachread,group=sample,label=1000*mean_perbp_eachread),vjust=0) +
      themebw + xlab("") + ylab("mutation/kb") + ggtitle("3.  Stats of [MIS/INS/DEL/MH] per kb of relevant isotypes\nmean(mean(each read's mut/each read's length in kb))")#\nSUM[i=0->i=n]{mutation[i]/bp[i]}/n => SUM( <Total mutation in a read> / <Total bp of that read> ) / <n = Total read)> => ( (1/70) + (4/65) + (6/60) + (2/90) ) / N ")
      if (length(grep("combine",dm1desc)) > 0) {
         p1.perbp_eachread = p1.perbp_eachread + geom_errorbar(aes(group=sample,ymin=1000*(mean_perbp_eachread-mean_perbp_eachread.se),ymax=1000*(mean_perbp_eachread+mean_perbp_eachread.se)),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=mean_perbp_eachread.pchar))
      }   
      
   pdf(paste("results1_",dm1desc,".pdf",sep=""),width=mysampletot*0.8,height=22)
   grid.arrange(p1tot,p1,p1.perbp_eachread,nrow=3,heights=c(1,4,4))
   dev.off()
   
   ## p2, breakdown mh by number
   p2.dm = dm1[grep("^(mis|mh|ins|del)$",dm1$type,invert=T,perl=T),]
   p2.dm = p2.dm[grep("(mis_[A-Z]..[A-Z]|(ins|del|mh)_[0-9])",p2.dm$type,perl=T),]

   p2tot = ggplot(p1.dm[p1.dm$type == "del",],aes(num,total_read)) +
      geom_bar(aes(fill=sample),color="black",stat="identity",position=position_dodge(width=0.9)) +
#      geom_errorbar(aes(group=sample,ymin=total_read-total_read.se,ymax=total_read+total_read.se),stat="identity",position=position_dodge(width=0.9),width=0.5) +
      facet_grid(myfill~.) +
      geom_text(aes(y=total_read,group=sample,label=total_read),stat="identity",position=position_dodge(width=0.9),vjust=0) +
      themebw + xlab("") + ylab("Number of read") + ggtitle("Total read for each sample and isotype")

   if (length(grep("combine",dm1desc)) > 0) {
      p2tot = p2tot + geom_errorbar(aes(group=sample,ymin=total_read-total_read.se,ymax=total_read+total_read.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=total_read.pchar))
   }   
   
   
   p2.dm.del = p2.dm[grep("del_[0-9]",p2.dm$type,perl=T),]
   #p2.dm.del = p2.dm.del[grep("mis_[A-Z]..[A-Z]$",p2.dm.del$type,perl=T),]
   p2.dm.del$num = gsub("^([0-9]+)\\..+$","\\1",p2.dm.del$sample,perl=T)
   p2.dm.del$num = as.numeric(as.character(p2.dm.del$num))
   p2.dm.del$x = gsub("^.+_([0-9]+)$","\\1",p2.dm.del$type,perl=T)
   p2.dm.del$x = as.numeric(as.character(p2.dm.del$x))
   
   p2.del = ggplot(p2.dm.del,aes(x,perc_atleast1)) +
      #geom_line(aes(color=sample)) +
      #geom_rect(aes(xmin=num-0.5,xmax=num+0.5,ymin=0,ymax=perc_atleast1,fill=sample),color="black") +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9),color="black") +
      facet_grid(myfill~.,scales = "free_y") +
      geom_text(aes(y=perc_atleast1,group=sample,label=perc_atleast1),vjust=0,stat="identity",position=position_dodge(width=0.9)) +
      themebwwithx + xlab("Length of Deletion (1, 2, 3, or 4 (4 = more than 3 bp))") + ylab("% of read") + ggtitle("% read with 1, 2, 3, or more than 3bp Deletion")
   if (length(grep("combine",dm1desc)) > 0) {
      p2.del = p2.del + geom_errorbar(aes(group=sample,ymin=perc_atleast1-perc_atleast1.se,ymax=perc_atleast1+perc_atleast1.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=perc_atleast1.pchar))
   }   
   
   
   p2.dm.ins = p2.dm[grep("ins_[0-9]",p2.dm$type,perl=T),]
   #p2.dm.ins = p2.dm.ins[grep("mis_[A-Z]..[A-Z]$",p2.dm.ins$type,perl=T),]
   #p2.dm.ins$myfill = paste(p2.dm.ins$treat,p2.dm.ins$isotype)#,"\nN =",p2.dm.ins$total_read)
   p2.dm.ins$num = gsub("^([0-9]+)\\..+$","\\1",p2.dm.ins$sample,perl=T)
   p2.dm.ins$num = as.numeric(as.character(p2.dm.ins$num))
   p2.dm.ins$x = gsub("^.+_([0-9]+)$","\\1",p2.dm.ins$type,perl=T)
   p2.dm.ins$x = as.numeric(as.character(p2.dm.ins$x))
   
   p2.ins = ggplot(p2.dm.ins,aes(x,perc_atleast1)) +
      #geom_line(aes(color=sample)) +
      #geom_rect(aes(xmin=num-0.5,xmax=num+0.5,ymin=0,ymax=perc_atleast1,fill=sample),color="black") +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9),color="black") +
      facet_grid(myfill~.,scales = "free_y") +
      geom_text(aes(y=perc_atleast1,group=sample,label=perc_atleast1),vjust=0,stat="identity",position=position_dodge(width=0.9)) +
      themebwwithx + xlab("Length of Insertion (1, 2, 3, or 4 (4 = more than 3 bp))") + ylab("% of read") + ggtitle("% read with 1, 2, 3, or more than 3bp Insertion")
   if (length(grep("combine",dm1desc)) > 0) {
      p2.ins = p2.ins + geom_errorbar(aes(group=sample,ymin=perc_atleast1-perc_atleast1.se,ymax=perc_atleast1+perc_atleast1.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=perc_atleast1.pchar))
   }   
   
   p2.dm.mh = p2.dm[grep("mh_[0-9]",p2.dm$type,perl=T),]
   #p2.dm.mh = p2.dm.mh[grep("mis_[A-Z]..[A-Z]$",p2.dm.mh$type,perl=T),]
   #p2.dm.mh$myfill = paste(p2.dm.mh$treat,p2.dm.mh$isotype)#,"\nN =",p2.dm.mh$total_read)
   p2.dm.mh$num = gsub("^([0-9]+)\\..+$","\\1",p2.dm.mh$sample,perl=T)
   p2.dm.mh$num = as.numeric(as.character(p2.dm.mh$num))
   p2.dm.mh$x = gsub("^.+_([0-9]+)$","\\1",p2.dm.mh$type,perl=T)
   p2.dm.mh$x = as.numeric(as.character(p2.dm.mh$x))
   
   p2.mh = ggplot(p2.dm.mh,aes(x,perc_atleast1)) +
      # geom_line(aes(color=sample)) +
      #geom_rect(aes(xmin=num-0.5,xmax=num+0.5,ymin=0,ymax=perc_atleast1,fill=sample),color="black") +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9),color="black") +
      facet_grid(myfill~.,scales = "free_y") +
      geom_text(aes(y=perc_atleast1,group=sample,label=perc_atleast1),vjust=0,stat="identity",position=position_dodge(width=0.9)) +
      themebwwithx + xlab("Length of MH (1, 2, 3, or 4 (4 = more than 3 bp))") + ylab("% of read") + ggtitle("% read with 1, 2, 3, or more than 3bp MH")
   if (length(grep("combine",dm1desc)) > 0) {
      p2.mh = p2.mh + geom_errorbar(aes(group=sample,ymin=perc_atleast1-perc_atleast1.se,ymax=perc_atleast1+perc_atleast1.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=perc_atleast1.pchar))
   }   
   
   p2.dm.mis = dm1[grep("mis_[A-Z]..[A-Z]$",dm1$type,perl=T),]
   #p2.dm.mis = p2.dm.mis[grep("mis_[A-Z]..[A-Z]$",p2.dm.mis$type,perl=T),]
   #p2.dm.mis$myfill = paste(p2.dm.mis$treat,p2.dm.mis$isotype)#,"\nN =",p2.dm.mis$total_read)
   p2.dm.mis$num = gsub("^([0-9]+)\\..+$","\\1",p2.dm.mis$sample,perl=T)
   p2.dm.mis$num = as.numeric(as.character(p2.dm.mis$num))
   
   p2.dm.mis$nuc1 = 1
   p2.dm.mis$nuc2 = 1
   p2.dm.mis[grep("^mis_G",p2.dm.mis$type,perl=T),]$nuc1 = 2
   p2.dm.mis[grep("^mis_C",p2.dm.mis$type,perl=T),]$nuc1 = 3
   p2.dm.mis[grep("^mis_T",p2.dm.mis$type,perl=T),]$nuc1 = 4
   p2.dm.mis[grep("^mis_.+G$",p2.dm.mis$type,perl=T),]$nuc2 = 2
   p2.dm.mis[grep("^mis_.+C$",p2.dm.mis$type,perl=T),]$nuc2 = 3
   p2.dm.mis[grep("^mis_.+T$",p2.dm.mis$type,perl=T),]$nuc2 = 4
   p2.dm.mis$order = as.numeric(as.character(paste(p2.dm.mis$nuc1,p2.dm.mis$nuc2,sep="")))
   p2.dm.mis$type = paste(p2.dm.mis$order,"\n",p2.dm.mis$type,sep="")
   
   p2.mis = ggplot(p2.dm.mis,aes(num,perc_atleast1)) +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9), color="black") +
      facet_grid(myfill~type,scales = "free_y") +
      geom_text(aes(y=perc_atleast1,group=sample,label=format(perc_atleast1,digits=2,scientific=F)),size=2,vjust=0) +
      themebwwithx + xlab("") + ylab("% of read with at least 1 [MIS_X->X]") + ggtitle("Stat of read [MIS_X->X] of relevant isotypes\nCOUNT[i=0->i=n]{mutation >= 1}/COUNT[i=0->i=n]{n} => <Total number of read with at least 1 [MIS_X->Y]>) / <Total read>")
   if (length(grep("combine",dm1desc)) > 0) {
      p2.mis = p2.mis + geom_errorbar(aes(group=sample,ymin=perc_atleast1-perc_atleast1.se,ymax=perc_atleast1+perc_atleast1.se),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=perc_atleast1.pchar))
   }   
   
   p2.mis.perbp_allread = ggplot(p2.dm.mis,aes(num,mean_perbp_allread)) +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9), color="black") +
      facet_grid(myfill~type,scales = "free_y") +
      geom_text(aes(y=mean_perbp_allread,group=sample,label=format(mean_perbp_allread,digits=2,scientific=F)),size=2,vjust=0) +
      themebwwithx + xlab("") + ylab("mutation/bp") + ggtitle("Stat of [MIS_X->X] relevant isotypes; mean(all reads)/mean(all bp)\nSUM[i=0->i=n]{mutation}/SUM[i=0->i=n]{bp} => <Total mutation of all reads> / <Total bp of all reads> => ( (1+4+6+2) / N ) / ( (70+65+60+90) / N )")
   if (length(grep("combine",dm1desc)) > 0) {
      p2.mis.perbp_allread = p2.mis.perbp_allread + geom_errorbar(aes(group=sample,ymin=mean_perbp_allread-mean_perbp_allread.se,ymax=mean_perbp_allread+mean_perbp_allread.se),stat="identity",position=position_dodge(width=0.9),width=0.5) + geom_text(aes(y=0,group=sample,label=mean_perbp_allread.pchar))
   }   
   
   p2.mis.perbp_eachread = ggplot(p2.dm.mis,aes(num,1000*mean_perbp_eachread)) +
      geom_bar(aes(fill=sample),stat="identity",position=position_dodge(width=0.9), color="black") +
      facet_grid(myfill~type,scales = "free_y") +
      geom_text(aes(y=1000*mean_perbp_eachread,group=sample,label=format(1000*mean_perbp_eachread,digits=2,scientific=F)),size=2,vjust=0) +
      themebw + xlab("") + ylab("mismatch mutation/kb") + ggtitle("Bar plot of mistmatch mutation/kb")
   #Stats of [MISMATCH X->Y] per bp of relevant isotypes; mean(mean(each read's mut/each read's bp))\nSUM[i=0->i=n]{mutation[i]/bp[i]}/n => SUM( <Total mutation in a read> / <Total bp of that read> ) / <n = Total read)> => ( (1/70) + (4/65) + (6/60) + (2/90) ) / N ")
   if (length(grep("combine",dm1desc)) > 0) {
      p2.mis.perbp_eachread = p2.mis.perbp_eachread + geom_errorbar(aes(group=sample,ymin=1000*(mean_perbp_eachread-mean_perbp_eachread.se),ymax=1000*(mean_perbp_eachread+mean_perbp_eachread.se)),stat="identity",position=position_dodge(width=0.9),width=0.5)+ geom_text(aes(y=0,group=sample,label=mean_perbp_eachread.pchar))
   }   
   
   
   pdf(paste("results2_",dm1desc,".pdf",sep=""),width=80, height=length(unique(as.character(dm1$myfill)))*3)
   grid.arrange(p2tot, p2.del, p2.ins, p2.mh, p2.mis.perbp_eachread,nrow=1,ncol=5,widths=c(1,1,1,1,4))
   dev.off()
}


mycombine = function(dm) {
   myse = function(x) {if(length(x) == 0) {return(0)} else {return(sd(x)/sqrt(length(x)))}}
   
   #total_read
   dm1 = as.data.frame(aggregate(dm$total_read,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))
   colnames(dm1) = c("sampleID","sample","treat","isotype","type","myfill","total_read")
   dm1$total_read.se = as.data.frame(aggregate(dm$total_read,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   head(dm1)
   
   #average_length
   dm1$average_length = as.data.frame(aggregate(dm$average_length,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$average_length.se = as.data.frame(aggregate(dm$average_length,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   
   #count_atleast1
   dm1$count_atleast1 = as.data.frame(aggregate(dm$count_atleast1,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$count_atleast1.se = as.data.frame(aggregate(dm$count_atleast1,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   
   #perc_atleast1
   dm1$perc_atleast1 = as.data.frame(aggregate(dm$perc_atleast1,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$perc_atleast1.se = as.data.frame(aggregate(dm$perc_atleast1,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   
   #count_total
   dm1$count_total = as.data.frame(aggregate(dm$count_total,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$count_total.se = as.data.frame(aggregate(dm$count_total,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   
   #mean_total
   dm1$mean_total = as.data.frame(aggregate(dm$mean_total,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$mean_total.se = as.data.frame(aggregate(dm$mean_total,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   
   #mean_perbp_allread
   dm1$mean_perbp_allread = as.data.frame(aggregate(dm$mean_perbp_allread,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$mean_perbp_allread.se = as.data.frame(aggregate(dm$mean_perbp_allread,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   
   #total_perbp_eachread
   dm1$total_perbp_allread = as.data.frame(aggregate(dm$total_perbp_eachread,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$total_perbp_eachread.se = as.data.frame(aggregate(dm$total_perbp_eachread,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   
   #mean_perbp_eachread
   dm1$mean_perbp_eachread = as.data.frame(aggregate(dm$mean_perbp_eachread,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),mean))$x
   dm1$mean_perbp_eachread.se = as.data.frame(aggregate(dm$mean_perbp_eachread,by=list(dm$sampleID, dm$sample,dm$treat,dm$isotype,dm$type,dm$myfill),myse))$x
   return(dm1)
}

find_p = function(dmtemp) {
   dmtemp$total_read.p = 1
   dmtemp$average_length.p = 1
   dmtemp$count_atleast1.p = 1
   dmtemp$perc_atleast1.p = 1
   dmtemp$count_total.p = 1
   dmtemp$mean_total.p = 1
   dmtemp$mean_perbp_allread.p = 1
   dmtemp$total_perbp_eachread.p = 1
   dmtemp$mean_perbp_eachread.p = 1
   dmtemp$total_read.pchar = ""
   dmtemp$average_length.pchar = ""
   dmtemp$count_atleast1.pchar = ""
   dmtemp$perc_atleast1.pchar = ""
   dmtemp$count_total.pchar = ""
   dmtemp$mean_total.pchar = ""
   dmtemp$mean_perbp_allread.pchar = ""
   dmtemp$total_perbp_eachread.pchar = ""
   dmtemp$mean_perbp_eachread.pchar = ""   
   typez = unique(dmtemp$type)
   myfillz = unique(dmtemp$myfill)
   samplez = unique(dmtemp$sample)
   for (i in 1:length(typez)) {
      print(paste(i,typez[i]))
      for (j in 1:length(myfillz)) {
         #         print(paste(i,j,typez[i],myfillz[j]))
         if (length(dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j],]$sample) > 0) {
            temp = dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j],]
            if (length(temp[grep("WT",temp$sample),]$sample) > 0) {
               
               WT = temp[grep("WT",temp$sample),]
               for (k in 1:length(samplez)) {
                  if (length(temp[grep(samplez[k],temp$sample),]$sample) > 0) {
                     TR = temp[grep(samplez[k],temp$sample),]
                     if (dim(WT)[1] > 1 | dim(TR)[1] > 1) {
                        myp = get_p(WT$total_read,TR$total_read)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$total_read.p = myp
                        }
                        myp = get_p(WT$average_length,TR$average_length)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$average_length.p = myp
                        }
                        myp = get_p(WT$count_atleast1,TR$count_atleast1)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$count_atleast1.p = myp
                        }
                        myp = get_p(WT$perc_atleast1,TR$perc_atleast1)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$perc_atleast1.p = myp
                        }
                        myp = get_p(WT$count_total,TR$count_total)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$count_total.p = myp
                        }
                        myp = get_p(WT$mean_total,TR$mean_total)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$mean_total.p = myp
                        }
                        myp = get_p(WT$mean_perbp_allread,TR$mean_perbp_allread)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$mean_perbp_allread.p = myp
                        }
                        myp = get_p(WT$total_perbp_eachread,TR$total_perbp_eachread)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$total_perbp_eachread.p = myp
                        }
                        myp = get_p(WT$mean_perbp_eachread,TR$mean_perbp_eachread)
                        if (!is.na(myp) & myp <= 0.05) {
                           dmtemp[dmtemp$type == typez[i] & dmtemp$myfill == myfillz[j] & dmtemp$sample == samplez[k],]$mean_perbp_eachread.p = myp
                        }
                     }
                  }
               }
            }
         }
      }
   }
   
   dmtemp$total_read.pchar = get_pchar(dmtemp$total_read.p)
   dmtemp$average_length.pchar = get_pchar(dmtemp$average_length.p)
   dmtemp$count_atleast1.pchar = get_pchar(dmtemp$count_atleast1.p)
   dmtemp$perc_atleast1.pchar = get_pchar(dmtemp$perc_atleast1.p)
   dmtemp$count_total.pchar = get_pchar(dmtemp$count_total.p)
   dmtemp$mean_total.pchar = get_pchar(dmtemp$mean_total.p)
   dmtemp$mean_perbp_allread.pchar = get_pchar(dmtemp$mean_perbp_allread.p)
   dmtemp$total_perbp_eachread.pchar = get_pchar(dmtemp$total_perbp_eachread.p)
   dmtemp$mean_perbp_eachread.pchar = get_pchar(dmtemp$mean_perbp_eachread.p)
   
   return(dmtemp)
}

get_p = function(WT,TR) {
   
   if (sd(WT) == 0 | sd(TR) == 0) {
      if (WT[1] == 0 & TR[1] == 0) {
         return(1)
      }
      if (WT[1] == 0) {
         add = sample(x = seq(0,TR[1]/1e1,by=TR[1]/1e2),size=length(WT),replace=FALSE)
         for (i in 1:length(WT)) {
            newWT = WT[i] + add[i]
#            print(paste(i,WT[i],"newWT = ",newWT))
            WT[i] = newWT
         }
      } else if (TR[1] == 0) {
         add = sample(x = seq(0,WT[1]/1e1,by=WT[1]/1e2),size=length(TR),replace=FALSE)
         for (i in 1:length(TR)) {
            newTR = TR[i] + add[i]
  #          print(paste(i,TR[i],"newTR = ",newTR))
            TR[i] = newTR
         }
      } else if (WT[1] < TR[1]) {
         add = sample(x = seq(0,WT[1]/1e1,by=WT[1]/1e2),size=length(WT),replace=FALSE)
         for (i in 1:length(WT)) {
            newWT = WT[i] + add[i]
  #         print(paste(i,WT[i],"newWT = ",newWT))
            WT[i] = newWT
         }
      } else if (TR[1] <= WT[1]) {
         add = sample(x = seq(0,TR[1]/1e1,by=TR[1]/1e2),size=length(TR),replace=FALSE)
         for (i in 1:length(TR)) {
            newTR = TR[i] + add[i]
#            print(paste(i,TR[i],"newTR = ",newTR))
            TR[i] = newTR
         }
      }
   }
         
   testtype = "greater"
   if (mean(WT) > mean(TR)) {
      testtype = "less"
   }
  # print(TR)
  # print(WT)
   myp = t.test(TR,WT,testtype)$p.value
   return(myp)
}

get_pchar = function(p) {
   pchar = rep("",length(p))
   if (length(p[p <= 0.05] > 0)) {
      pchar[p <= 0.05] = "*"
   }
   if (length(p[p <= 0.01] > 0)) {
      pchar[p <= 0.01] = "**"
   }
   if (length(p[p <= 0.001] > 0)) {
      pchar[p <= 0.001] = "***"
   }
   return(pchar)   
}


file1 = "/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/results.csv"
file2 = "/group/stella/Work/Data/Fastq/200903_SLIMS4168_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/results.csv"
dm1 = myparse(file1,"run1_SLIMS0323")
dm2 = myparse(file2,"run2_SLIMS4168")
dmbackup = mycombine(rbind(dm1,dm2))
dmtempbackup = find_p(rbind(dm1,dm2))
dmtemp = dmtempbackup[grep("run1_",dmtempbackup$file),]
dmtempz = subset(dmtemp,select=c(-total_read,-average_length,-count_atleast1,-perc_atleast1,-count_total,-mean_total,-mean_perbp_allread,-total_perbp_eachread,-mean_perbp_eachread))
dm = merge(dmbackup,dmtempz,by=c("sampleID","sample","treat","isotype","type","myfill"),all=T)
# file1

dm1.IgG1and3 = dm1[grep("(IgG1.+IgG1|IgG3.+IgG3|IgG1.+IgG3|01. IgM)",paste(dm1$isotype,dm1$treat),perl=T),]
if (dim(dm1.IgG1and3)[1] > 0) {
   mygraph(dm1.IgG1and3,"run1_IgG1_and_IgG3_SLIMS0323")
}

# file2
dm2.IgG1and3 = dm2[grep("(IgG1.+IgG1|IgG3.+IgG3|IgG1.+IgG3|01. IgM)",paste(dm2$isotype,dm2$treat),perl=T),]
if (dim(dm2.IgG1and3)[1] > 0) {
   mygraph(dm2.IgG1and3,"run2_IgG1_and_IgG3_SLIMS4168")
}

# combined
dm.IgG1and3 = dm[grep("(IgG1.+IgG1|IgG3.+IgG3|IgG1.+IgG3|01. IgM)",paste(dm$isotype,dm$treat),perl=T),]
if (dim(dm.IgG1and3)[1] > 0) {
   mygraph(dm.IgG1and3,"run1_IgG1_and_IgG3_combined")
}

# big picture
mygraph(dm1,"run1_ALL_SLIMS0323")
mygraph(dm2,"run2_ALL_SLIM4168")
mygraph(dm,"run1_ALL_combined")


parse_juncdist = function(x) {
   files1 = dir("/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/",".juncdist$")
   files1desc = "run1_SLIMS0323"
   files2 = dir("/group/stella/Work/Data/Fastq/200903_SLIMS4168_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/",".juncdist$")
   files2desc = "run2_SLIMS4168"
   
   for (i in 1:length(files1)) {
      file1 = paste("/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/",files1[i],sep="")
      print(paste("1.",i,file1))
      if (i == 1) {
         juncdm1 = myparse(file1, files1desc)
      } else {
         juncdm1 = rbind(juncdm1, myparse(file1, files1desc))
      }
   }
   #colnames(juncdm1)[colnames(juncdm1) == "strand2"] = "strand"
   juncdm1$one = 1

   # juncdm1total = aggregate(juncdm1$one,by=list(juncdm1$sampleID, juncdm1$sample, juncdm1$treat, juncdm1$isotype, juncdm1$strand, juncdm1$myfill),sum)
   # colnames(juncdm1total) = c("sampleID","sample","treat","isotype","strand","myfill","total_read")
   # juncdm1 = merge(juncdm1,juncdm1total,by=c("sampleID","sample","treat","isotype","strand","myfill"))
   juncdm1total = aggregate(juncdm1$one,by=list(juncdm1$sampleID, juncdm1$sample, juncdm1$treat, juncdm1$isotype,  juncdm1$myfill),sum)
   colnames(juncdm1total) = c("sampleID","sample","treat","isotype","myfill","total_read")
   juncdm1total$total_read =juncdm1total$total_read /2 
   juncdm1 = merge(juncdm1,juncdm1total,by=c("sampleID","sample","treat","isotype","myfill"))
   juncdm1$juncperc = as.integer(juncdm1$juncperc)
   if (length(juncdm1[juncdm1$juncperc > 100,]$juncperc) > 0) {
      juncdm1[juncdm1$juncperc > 100,]$juncperc = 100
   }
   if (length(juncdm1[juncdm1$juncdist > 15000,]$juncdist > 0)) {
      print(paste("juncdm1, dist1 total =",length(juncdm1[juncdm1$junc2dist > 15000,]$junc2dist)))
      juncdm1[juncdm1$juncdist > 15000,]$juncdist = 15000 
   }
   juncdm1$juncdist = as.integer(juncdm1$juncdist / 150)

   for (i in 1:length(files2)) {
      file2 = paste("/group/stella/Work/Data/Fastq/200903_SLIMS4168_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/",files2[i],sep="")
      print(paste("2.",i,file2))
      if (i == 1) {
         juncdm2 = myparse(file2, files2desc)
      } else {
         juncdm2 = rbind(juncdm2, myparse(file2, files2desc))
      }
   }
   colnames(juncdm2)[colnames(juncdm2) == "strand2"] = "strand"
   juncdm2$one = 1
   # juncdm2total = aggregate(juncdm2$one,by=list(juncdm2$sampleID, juncdm2$sample, juncdm2$treat, juncdm2$isotype, juncdm2$strand, juncdm2$myfill),sum)
   # colnames(juncdm2total) = c("sampleID","sample","treat","isotype","strand","myfill","total_read")
   # juncdm2 = merge(juncdm2,juncdm2total,by=c("sampleID","sample","treat","isotype","strand","myfill"))
   juncdm2total = aggregate(juncdm2$one,by=list(juncdm2$sampleID, juncdm2$sample, juncdm2$treat, juncdm2$isotype, juncdm2$myfill),sum)
   colnames(juncdm2total) = c("sampleID","sample","treat","isotype","myfill","total_read")
   juncdm2total$total_read =juncdm2total$total_read /2
   juncdm2 = merge(juncdm2,juncdm2total,by=c("sampleID","sample","treat","isotype","myfill"))
   juncdm2$juncperc = as.integer(juncdm2$juncperc)
   if (length(juncdm2[juncdm2$juncperc > 100,]$juncperc) > 0) {
      juncdm2[juncdm2$juncperc > 100,]$juncperc = 100
   }
   if (length(juncdm2[juncdm2$juncdist > 15000,]$juncdist > 0)) {
      print(paste("juncdm2, dist1 total =",length(juncdm2[juncdm2$juncdist > 15000,]$juncdist)))
      juncdm2[juncdm2$juncdist > 15000,]$juncdist = 15000 
   }
   juncdm2$juncdist = as.integer(juncdm2$juncdist / 150)
   
   juncdm1backup = juncdm1
   juncdm2backup = juncdm2
   
   #by perc len
   juncdm1 = juncdm1backup
   juncdm1temp1 = aggregate(juncdm1$one,by=list(juncdm1$sampleID, juncdm1$sample, juncdm1$treat, juncdm1$isotype, juncdm1$type, juncdm1$strand, juncdm1$myfill,juncdm1$total_read,juncdm1$juncperc,juncdm1$junclen),sum)
   colnames(juncdm1temp1) = c("sampleID", "sample", "treat", "isotype", "type","strand","myfill","total_read","xpos","junclen","count")
   juncdm1 = juncdm1temp1
   juncdm1$perc = as.integer(juncdm1$count / juncdm1$total_read * 10000+0.5)/100
   juncdm1$filedesc = files1desc
   juncdm1$category = "by_perc_len"
   juncdm1.byperclen = juncdm1
   
   juncdm2 = juncdm2backup
   juncdm2temp1 = aggregate(juncdm2$one,by=list(juncdm2$sampleID, juncdm2$sample, juncdm2$treat, juncdm2$isotype, juncdm2$type, juncdm2$strand, juncdm2$myfill,juncdm2$total_read,juncdm2$juncperc,juncdm2$junclen),sum)
   colnames(juncdm2temp1) = c("sampleID", "sample", "treat", "isotype", "type","strand","myfill","total_read","xpos","junclen","count")
   juncdm2 = juncdm2temp1
   juncdm2$perc = as.integer(juncdm2$count / juncdm2$total_read * 10000+0.5)/100
   juncdm2$filedesc = files2desc
   juncdm2$category = "by_perc_len"
   juncdm2.byperclen = juncdm2
   
   #by 50bp bin
   juncdm1 = juncdm1backup
   juncdm1temp1 = aggregate(juncdm1$one,by=list(juncdm1$sampleID, juncdm1$sample, juncdm1$treat, juncdm1$isotype, juncdm1$type, juncdm1$strand, juncdm1$myfill,juncdm1$total_read,juncdm1$juncdist,juncdm1$junclen),sum)
   colnames(juncdm1temp1) = c("sampleID", "sample", "treat", "isotype", "type","strand","myfill","total_read","xpos","junclen","count")
   juncdm1 = juncdm1temp1
   juncdm1$perc = as.integer(juncdm1$count / juncdm1$total_read * 10000+0.5)/100
   juncdm1$filedesc = files1desc
   juncdm1$category = "by_150bp_bin"
   juncdm1.by50bpbin = juncdm1
   
   juncdm2 = juncdm2backup
   juncdm2temp1 = aggregate(juncdm2$one,by=list(juncdm2$sampleID, juncdm2$sample, juncdm2$treat, juncdm2$isotype,juncdm2$type,  juncdm2$strand, juncdm2$myfill,juncdm2$total_read,juncdm2$juncdist,juncdm2$junclen),sum)
   colnames(juncdm2temp1) = c("sampleID", "sample", "treat", "isotype", "type","strand","myfill","total_read","xpos","junclen","count")
   juncdm2 = juncdm2temp1
   juncdm2$perc = as.integer(juncdm2$count / juncdm2$total_read * 10000+0.5)/100
   juncdm2$filedesc = files2desc
   juncdm2$category = "by_150bp_bin"
   juncdm2.by50bpbin = juncdm2
   
   
   juncdm = rbind(juncdm1.byperclen, juncdm2.byperclen, juncdm1.by50bpbin, juncdm2.by50bpbin)
   juncdm1total$filedesc = files1desc
   juncdm2total$filedesc = files2desc
   juncdmtotal = rbind(juncdm1total,juncdm2total)
   return(juncdm)
}

myresbackup = parse_juncdist()
myresbackup[myresbackup$strand == "-",]$perc = myresbackup[myresbackup$strand == "-",]$perc * -1 

myres = myresbackup

categories = unique(myres$category)
filedesc = unique(myres$filedesc)
typez = unique(myres$type)
temp = myres[grep("(IgG1.+IgG1$|IgG1.+IgM$|IgG3.+IgG1$|IgG3.+IgG3$|IgG3.+IgM$)",paste(myres$treat,myres$isotype),perl=T),]
temp = subset(temp,select=-sampleID)
temp$myfillz = paste(temp$category, temp$filedesc, temp$sample, temp$treat, temp$isotype, temp$strand, temp$myfill, temp$total_read, temp$junclen, temp$type)
myfillz = unique(temp$myfillz)
for (i in 1:length(categories)) {
   if (categories[i] == "by_perc_len") {
      xposall = seq(0,100,by=1)
   } else if (categories[i] == "by_150bp_bin") {
      xposall = seq(0,15000/150,by=150/150)
   }
   for (k in 1:length(myfillz)) {
      if (k %% 10 == 0) {
         print(paste(categories[i],i,k))
      }
      tempall = temp[temp$myfillz == myfillz[k] & temp$category == categories[i],]
      a = tempall[1,]
      tempxpos = xposall[!xposall %in% unique(tempall$xpos)]
      if (length(tempxpos) > 0) {
         temp2 = data.frame(category = a$category, filedesc=a$filedesc, sample=a$sample, treat=a$treat, isotype=a$isotype, strand=a$strand, myfill=a$myfill, total_read=a$total_read, junclen=a$junclen, type=a$type,xpos=tempxpos,perc = 0,myfillz = myfillz[k],count=0)
         temp = rbind(temp,temp2)
      }
   }
}

a = unique(temp[temp$category != "by_perc_len",]$xpos)
a = a[order(a)]
a
#      temp =merge(temp,temp3,by=c("myfillz","xpos"),all=T)
myres = temp
myres = myres[grep("(IgG1.+IgG1$|IgG1.+IgM$|IgG3.+IgG1$|IgG3.+IgG3$|IgG3.+IgM$)",paste(myres$treat,myres$isotype),perl=T),]

categories = unique(myres$category)
filedesc = unique(myres$filedesc)
typez = unique(myres$type)


dm = myres
dm[dm$perc < 0,]$perc =dm[dm$perc < 0,]$perc *-1 
dm = dm[dm$category == "by_perc_len",]
dm$cat = "Last 36%"
dm[dm$xpos > (100*36/100),]$cat = "First 63%"
dm = dm[dm$cat == "Last 36%",]
dms = aggregate(dm$perc,by=list(dm$myfill,dm$category,dm$filedesc,dm$type,dm$sample,dm$cat),sum)
colnames(dms) = c("myfill","category","filedesc","type","sample","xpos","perc")

pdf("results3_barplot.pdf",width=10,height=20)
for (i in 1:length(typez)) {
   temp = dms[dms$type== typez[i],]   
   p1 = ggplot(temp,aes(xpos,perc)) +
      geom_bar(aes(fill=sample),stat="identity",position="dodge",color="black") +
      geom_text(aes(group=sample,label=perc),stat="identity",position=position_dodge(width=0.9)) +
      theme_bw() + theme(panel.grid=element_blank()) +
      facet_grid(myfill~filedesc) +
      ylab("% of total read") + xlab(paste("Last 36% of length between [exon1 ATG-5000bp] until [last exon] of Ig gene",typez[i]))
   print(p1)
}
dev.off()

pdf("results3.pdf",height=8,width=10)
for (i in 1:length(categories)) {
   for (j in 1:length(filedesc)) {
      for (k in 1:length(typez)) {
         if (categories[i] == "by_perc_len") {
            temp = myres[myres$category == categories[i] & myres$filedesc == filedesc[j] & myres$type == typez[k],]
            print(paste(i,j,k,categories[i],filedesc[j],typez[k],dim(temp)[1]))
            if(length(temp$perc) > 0) {
               temptotal = temp[temp$xpos == 0 & temp$strand == "-",]
               temptotal$y = 2
               p1 = ggplot(temp,aes(xpos,perc)) +
                  geom_line(aes(color=strand),size=0.2) +
                  geom_text(data=temptotal,aes(x=1,y=y,label=paste("n = ",total_read,sep="")),hjust=0) +
                  facet_grid(myfill~sample) +
                  theme_bw() + theme(panel.grid = element_blank()) +
                  ggtitle(paste(categories[i],filedesc[j],typez[k])) + ylab("% of total read") + xlab("% length from -5kb of ATG to 3'end of gene") +
                  scale_color_manual(values=c("+"="red2","-"="blue2"))
               if (categories[i] == "by_perc_len") {
                  if (typez[k] == "BAIT") {
                     p1 = p1 + coord_cartesian(xlim=c(0,100),ylim=c(-25,25))
                     p1 = p1 + annotate(geom="segment",x=25,xend=25,y=-25,yend=25,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=50,xend=50,y=-25,yend=25,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=75,xend=75,y=-25,yend=25,lty=2,alpha=0.25)
                  } else {
                     p1 = p1 + coord_cartesian(xlim=c(0,100),ylim=c(-10,10))
                     p1 = p1 + annotate(geom="segment",x=25,xend=25,y=-10,yend=10,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=50,xend=50,y=-10,yend=10,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=75,xend=75,y=-10,yend=10,lty=2,alpha=0.25)
                  }
               } else {
                  if (typez[k] == "BAIT") {
                     p1 = p1 + coord_cartesian(xlim=c(0,15000/150),ylim=c(-25,25))
                     p1 = p1 + annotate(geom="segment",x=25,xend=25,y=-25,yend=25,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=50,xend=50,y=-25,yend=25,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=75,xend=75,y=-25,yend=25,lty=2,alpha=0.25)
                  } else {
                     p1 = p1 + coord_cartesian(xlim=c(0,15000/150),ylim=c(-10,10))
                     p1 = p1 + annotate(geom="segment",x=25,xend=25,y=-10,yend=10,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=50,xend=50,y=-10,yend=10,lty=2,alpha=0.25)
                     p1 = p1 + annotate(geom="segment",x=75,xend=75,y=-10,yend=10,lty=2,alpha=0.25)
                  }
               }
               
               print(p1)
               # p2 = ggplot(temptotal,aes(xpos,total_read)) +
               #    geom_bar(aes(fill=strand),stat="identity",position="dodge") +
               #    facet_grid(myfill~sample) +
               #    theme_bw() + theme(panel.grid = element_blank()) +
               #    ggtitle(paste(categories[i],filedesc[j],typez[k])) + ylab("% of total read") + xlab("") +
               #    scale_fill_manual(values=c("+"="red2","-"="blue2"))
            }
         }
      }
   }
}
dev.off()

