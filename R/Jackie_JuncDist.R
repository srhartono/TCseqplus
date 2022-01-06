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
myres = myres[grep("(IgG1.+IgG1$|IgG3.+IgG3$)",paste(myres$treat,myres$isotype),perl=T),]

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

pdf("results3_barplot.pdf",width=18,height=7)
for (i in 1:length(typez)) {
   temp = dms[dms$type== typez[i],]   
   temp = temp[grep("IgM",temp$myfill,invert = T),]
   temp$myfill2 = temp$myfill
   temp$myfill2 = gsub("^(.+)undig(.+)$","A. NOT CUTMUHSPACE\\1\\2",temp$myfill2,perl=T)
   temp$myfill2 = gsub("^(.+)diges(.+)$","B. CUTMUHSPACE\\1\\2",temp$myfill2,perl=T)
   temp$myfill2 = gsub("^(.+)([0-9]+)\\. (Ig[A-Z]+[0-9]+)  (germ|rich) (0[0-9]\\..+)$","\\1GROUP \\2 (\\3) \\4MUHSPACEPREY at \\5",temp$myfill2,perl=T)
   temp$myfill3 = temp$myfill2
   temp$myfill2 = gsub(" germ","",temp$myfill2,perl=T)
   temp$myfill2 = gsub(" rich","",temp$myfill2,perl=T)
   temp$myfill3 = gsub(" germ","",temp$myfill3,perl=T)
   temp$myfill3 = gsub(" rich","",temp$myfill3,perl=T)
   temp$myfill2 = gsub("GROUP 2","GROUP 1",temp$myfill2,perl=T)
   temp$myfill4 = gsub("(A|B)\\. (NOT CUT|CUT)","",temp$myfill2,perl=T)
   temp$myfill4 = gsub("GROUP 2","GROUP 1",temp$myfill4,perl=T)
   temp$myfill2 = gsub("MUHSPACE","\\\n",temp$myfill2,perl=T)
   temp$myfill3 = gsub("MUHSPACE","\\\n",temp$myfill3,perl=T)
   temp$myfill4 = gsub("MUHSPACE","\\\n",temp$myfill4,perl=T)
   
   temp2 = as.data.frame(aggregate(temp$perc,by=list(temp$category,temp$type,temp$sample,temp$xpos,temp$myfill2),mean))
   colnames(temp2) = c("category","type","sample","xpos","myfill2","mymean")
   temp2$myse = as.data.frame(aggregate(temp$perc,by=list(temp$category,temp$type,temp$sample,temp$xpos,temp$myfill2),function(x)sd(x)/sqrt(length(x))))$x
   temp2$mymean = as.integer(temp2$mymean * 100)/100
   temp2$myymax = as.integer((temp2$mymean + temp2$myse)*100)/100
   temp2$myymin = as.integer((temp2$mymean - temp2$myse)*100)/100

   temp4 = as.data.frame(aggregate(temp$perc,by=list(temp$category,temp$type,temp$sample,temp$xpos,temp$myfill4),mean))
   colnames(temp4) = c("category","type","sample","xpos","myfill4","mymean")
   temp4$myse = as.data.frame(aggregate(temp$perc,by=list(temp$category,temp$type,temp$sample,temp$xpos,temp$myfill4),function(x)sd(x)/sqrt(length(x))))$x
   temp4$mymean = as.integer(temp4$mymean * 100)/100
   temp4$myymax = as.integer((temp4$mymean + temp4$myse)*100)/100
   temp4$myymin = as.integer((temp4$mymean - temp4$myse)*100)/100

   myp2 = data.frame()
   myfill2z = unique(temp2$myfill2)
   samplez = unique(temp2$sample)
   for (j in 1:length(myfill2z)) {
  #    print(paste("Doing",j,myfill2z[j]))
      tempmyp2 = temp2[temp2$myfill2 == myfill2z[j] & temp2$sample == "1. WT",]
      if (dim(tempmyp2)[1] > 0) {
         print(head(tempmyp2))
         tempmyp2$p = 1
         tempmyp2$pchar = ""
         if (j == 1) {
            myp2 = tempmyp2
         } else {
            myp2 = rbind(myp2, tempmyp2)
         }
      }
      for (k in 2:length(samplez)) {
         print(paste("myp2",j,myfill2z[j],k,samplez[k]))
         WT = temp[temp$myfill2 == myfill2z[j] & temp$sample == "1. WT",]
         TR = temp[temp$myfill2 == myfill2z[j] & temp$sample == samplez[k],]
         if (dim(WT)[1] > 0 & dim(TR)[1] > 0) {
            WTperc = c()
            TRperc = c()
            if (dim(WT)[1] > 0) {
               WTperc = WT$perc
            }
            if (dim(TR)[1] > 0) {
               TRperc = TR$perc
            }
            tempp = get_p(WTperc,TRperc)
            temppchar = get_pchar(tempp)
            tempmyp2 = temp2[temp2$myfill2 == myfill2z[j] & temp2$sample == samplez[k],]
            tempmyp2$p = tempp
            tempmyp2$pchar = temppchar
            print(head(tempmyp2))
            myp2 = rbind(myp2, tempmyp2)
         }
      }
   }
   
   print("myp4")
   myp4 = data.frame()
   myfill4z = unique(temp4$myfill4)
   samplez = unique(temp4$sample)
   for (j in 1:length(myfill4z)) {
      tempmyp4 = temp4[temp4$myfill4 == myfill4z[j] & temp4$sample == "1. WT",]
      if (dim(tempmyp4)[1] > 0) {
         tempmyp4$p = 1
         tempmyp4$pchar = ""
         if (j == 1) {
            myp4 = tempmyp4
         } else {
            myp4 = rbind(myp4, tempmyp4)
         }
      }
      for (k in 2:length(samplez)) {
         print(paste("myp4",j,myfill4z[j],k,samplez[k]))
         WT = temp[temp$myfill4 == myfill4z[j] & temp$sample == "1. WT",]
         TR = temp[temp$myfill4 == myfill4z[j] & temp$sample == samplez[k],]
         if (dim(WT)[1] > 0 & dim(TR)[1] > 0) {
            WTperc = c()
            TRperc = c()
            if (dim(WT)[1] > 0) {
               WTperc = WT$perc
            }
            if (dim(TR)[1] > 0) {
               TRperc = TR$perc
            }
            tempp = get_p(WTperc,TRperc)
            temppchar = get_pchar(tempp)
            tempmyp4 = temp4[temp4$myfill4 == myfill4z[j] & temp4$sample == samplez[k],]
            tempmyp4$p = tempp
            tempmyp4$pchar = temppchar
            
            myp4 = rbind(myp4, tempmyp4)
         }
      }
   }
   
      #p1 = ggplot(temp2,aes(xpos,perc)) +
   p1 = ggplot(temp4,aes(xpos,mymean)) +
      geom_bar(aes(fill=sample),stat="identity",position="dodge",color="black") +
      geom_errorbar(aes(group=sample,ymin=myymin,ymax=myymax),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
      geom_text(aes(group=sample,y=mymean+0.5,label=mymean),stat="identity",position=position_dodge(width=0.9)) +
      geom_text(data=myp4,aes(group=sample,y=0,label=pchar),stat="identity",position=position_dodge(width=0.9),size=5) +
      geom_text(data=myp4,aes(group=sample,y=-0.5,label=paste("p=",as.integer(p*1000)/1000,sep="")),stat="identity",position=position_dodge(width=0.9),size=2) +
      theme_bw() + theme(panel.grid=element_blank()) +
      #facet_grid(myfill~filedesc) +
      facet_grid(.~myfill4) +
      ylab("% of total read") + xlab(paste("Last 36% of length between [exon1 ATG-5000bp] until [last exon] of Ig gene",typez[i]))
   print(p1)
   
   p2 = ggplot(temp2,aes(xpos,mymean)) +
      geom_bar(aes(fill=sample),stat="identity",position="dodge",color="black") +
      geom_errorbar(aes(group=sample,ymin=myymin,ymax=myymax),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
      geom_text(aes(group=sample,y=mymean+0.5,label=mymean),stat="identity",position=position_dodge(width=0.9),) +
      geom_text(data=myp2,aes(group=sample,y=0,label=pchar),stat="identity",position=position_dodge(width=0.9),size=5) +
      geom_text(data=myp2,aes(group=sample,y=-0.5,label=paste("p=",as.integer(p*1000)/1000,sep="")),stat="identity",position=position_dodge(width=0.9),size=2) +
      theme_bw() + theme(panel.grid=element_blank()) +
      #facet_grid(myfill~filedesc) +
      facet_grid(.~myfill2) +
      ylab("% of total read") + xlab(paste("Last 36% of length between [exon1 ATG-5000bp] until [last exon] of Ig gene",typez[i]))
   print(p2)
   
   p3 = ggplot(temp,aes(xpos,perc)) +
      geom_bar(aes(fill=sample),stat="identity",position="dodge",color="black") +
      #geom_errorbar(aes(group=sample,ymin=myymin,ymax=myymax),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
      geom_text(aes(group=sample,y=perc+0.5,label=perc),stat="identity",position=position_dodge(width=0.9)) +
      theme_bw() + theme(panel.grid=element_blank()) +
      #facet_grid(myfill~filedesc) +
      facet_grid(filedesc~myfill3) +
      ylab("% of total read") + xlab(paste("Last 36% of length between [exon1 ATG-5000bp] until [last exon] of Ig gene",typez[i]))
   print(p3)
}
dev.off()

pdf("results3.pdf",height=10,width=20)
for (i in 1:length(categories)) {
   for (j in 1:length(filedesc)) {
      for (k in 1:length(typez)) {
         if (categories[i] == "by_perc_len") {
            temp = myres[myres$category == categories[i] & myres$filedesc == filedesc[j] & myres$type == typez[k],]
            
            
            temp$myfill2 = temp$myfill
            temp$myfill2 = gsub("^(.+)undig(.+)$","A. NOT CUTMUHSPACE\\1\\2",temp$myfill2,perl=T)
            temp$myfill2 = gsub("^(.+)diges(.+)$","B. CUTMUHSPACE\\1\\2",temp$myfill2,perl=T)
            temp$myfill2 = gsub("^(.+)([0-9]+)\\. (Ig[A-Z]+[0-9]+)  (germ|rich) (0[0-9]\\..+)$","\\1GROUP \\2 (\\3) \\4MUHSPACEPREY at \\5",temp$myfill2,perl=T)
            temp$myfill3 = temp$myfill2
            temp$myfill2 = gsub(" germ","",temp$myfill2,perl=T)
            temp$myfill2 = gsub(" rich","",temp$myfill2,perl=T)
            temp$myfill3 = gsub(" germ","",temp$myfill3,perl=T)
            temp$myfill3 = gsub(" rich","",temp$myfill3,perl=T)
            temp$myfill2 = gsub("GROUP 2","GROUP 1",temp$myfill2,perl=T)
            temp$myfill4 = gsub("(A|B)\\. (NOT CUT|CUT)","",temp$myfill2,perl=T)
            temp$myfill4 = gsub("GROUP 2","GROUP 1",temp$myfill4,perl=T)
            temp$myfill2 = gsub("MUHSPACE","\\\n",temp$myfill2,perl=T)
            temp$myfill3 = gsub("MUHSPACE","\\\n",temp$myfill3,perl=T)
            temp$myfill4 = gsub("MUHSPACE","\\\n",temp$myfill4,perl=T)
            
            print(paste(i,j,k,categories[i],filedesc[j],typez[k],dim(temp)[1]))
            if(length(temp$perc) > 0) {
               temptotal = temp[temp$xpos == 0 & temp$strand == "-",]
               temptotal$y = 2
               p1 = ggplot(temp,aes(xpos,perc)) +
                  annotate(geom="segment",x=0,xend=100,y=0,yend=0,alpha=0.25) +
                  geom_line(aes(color=strand),size=0.2) +
                  facet_grid(sample~myfill3) +
                  theme_bw() + theme(panel.grid = element_blank()) +
                  ggtitle(paste(categories[i],filedesc[j],typez[k])) + ylab("% of total read") + xlab("% length from -5kb of ATG to 3'end of gene\n0 = 3'end of gene\n100 = -5kb of ATG\nRed box = Last 36% of gene length, to detect increased distal junctions)") +
                  scale_color_manual(values=c("+"="red2","-"="blue2"))
               
               if (typez[k] == "BAIT") {
                  p1 = p1 + geom_text(data=temptotal,aes(x=80,y=16,label=paste("n = ",total_read,sep="")),hjust=0)
                  p1 = p1 + annotate(geom="rect",xmin=0,xmax=36,ymin=-20,ymax=20,lty=2,alpha=0.25,color="red2",fill=rgb(1,1,1,alpha = 0))
               } else {
                  p1 = p1 + geom_text(data=temptotal,aes(x=80,y=8,label=paste("n = ",total_read,sep="")),hjust=0)
                  p1 = p1 + annotate(geom="rect",xmin=0,xmax=36,ymin=-10,ymax=10,lty=2,alpha=0.25,color="red2",fill=rgb(1,1,1,alpha = 0))
               }
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

