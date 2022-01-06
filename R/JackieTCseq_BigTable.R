library(ggplot2)
library(reshape2)
library(grid)
library(gridExtra)

se = function(x){sd(x)/sqrt(length(x))}

myparse = function(file,filedesc) {
   dfall = read.table(file,sep=",",header=T)[,-1]
   #   dfall = subset(dfall,select=-myfill2)
   for (i in 1:8) {
      dfall[,i] = as.character(dfall[,i])
   }
   dfall$mypaste = dfall[,1]
   for (i in 1:8) {
      dfall$mypaste = paste(dfall$mypaste,dfall[,i])
      dfall$mypaste = as.character(dfall$mypaste)
   }
   dfall2 = dfall[,seq(1,8)]
   dfall2$mypaste = dfall2[,1]
   for (i in 1:8) {
      dfall2$mypaste = paste(dfall2$mypaste,dfall2[,i])
      dfall2$mypaste = as.character(dfall2$mypaste)
   }
   dfall = dfall[,seq(9,dim(dfall)[2])]
   dfall = cbind(dfall2$mypaste,dfall)
   dfall = dfall[,-dim(dfall)[2]]
   colnames(dfall)[1] = "mypaste"
   # dfall$sampleID = as.character(dfall$sampleID)
   # dfall$sample = as.character(dfall$sample)
   # dfall$treat = as.character(dfall$treat)
   # dfall$isotype = as.character(dfall$isotype)
   mypastez = unique(dfall2$mypaste)
   for (i in 1:length(mypastez)) {
      mypaste = mypastez[i]
      mymean = apply(dfall[dfall$mypaste == mypastez[i],-1],2,mean)
      mytotal = dim(dfall[dfall$mypaste == mypastez[i],-1])[1]
      myse = apply(dfall[dfall$mypaste == mypastez[i],-1],2,se)
      colnamez = colnames(dfall)[-1]
      dm1temp = data.frame(mymean=mymean,myse=myse)
      for (j in 1:dim(dfall2)[2]) {
         dm1temp[,(dim(dm1temp)[2]+1)] = unique(dfall2[dfall2$mypaste == mypastez[i],j])[1]
         colnames(dm1temp)[dim(dm1temp)[2]] = colnames(dfall2)[j]
      }
      dm1temp$cat = rownames(dm1temp)
      mytotaltemp = dm1temp[1,]
      mytotaltemp$cat = "total_read"
#      print(paste(i,mypastez[i],"total read = ",mytotal))
      mytotaltemp$mymean = mytotal
      mytotaltemp$myse = 0
      dm1temp = rbind(dm1temp,mytotaltemp)
      rownames(dm1temp) = seq(1,dim(dm1temp)[1])
      if (i == 1) {
         dm1 = dm1temp
      } else {
         dm1 = rbind(dm1,dm1temp)
      }
   }
   
   dm1$filedesc = filedesc
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
   if (is.na(WT[1])) {
      WT[1] = 0
   }
   WTsd = sd(WT)
   if (is.na(WTsd)) {
      WTsd = 0
   }
   
   TRsd = sd(TR)
   if(is.na(TRsd)) {
      TRsd = 0
   }
   if (WTsd + TRsd == 0) {
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
   if (myp > 0.5) {
      myp = 1 - myp
   }
   return(myp)
}

get_pchar = function(p) {
   pchar = rep("",length(p))
   if (length(is.na(p)) > 0) {
      pchar[is.na(p)] = ""
   }
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



file1 = paste("/group/stella/Work/Data/Fastq/200323_SLIMS0323_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/both.csv",sep="")
file2 = paste("/group/stella/Work/Data/Fastq/200903_SLIMS4168_BC_mm_TCseq_JackieBarlowTCseq/2_Peak/results/both.csv",sep="")

dm1 = as.data.frame(myparse(file1,"rep1"))
dm2 = as.data.frame(myparse(file2,"rep2"))

dm1backup = dm1
dm2backup = dm2

dm1[grep("(^len$|total_read|_perkb_)",dm1$cat,invert=T,perl=T),]$mymean = dm1[grep("(^len$|total_read|_perkb_)",dm1$cat,invert=T,perl=T),]$mymean * 100
dm1[grep("(^len$|total_read|_perkb_)",dm1$cat,invert=T,perl=T),]$myse = dm1[grep("(^len$|total_read|_perkb_)",dm1$cat,invert=T,perl=T),]$myse * 100
# dm1[grep("\\\\nCUT\\\\n",dm1$myfill2,perl=T),]$myfill2 = gsub("CUT","2 CUT",dm1[grep("\\\\nCUT\\\\n",dm1$myfill2,perl=T),]$myfill2,perl=T)
# dm1[grep("\\\\nUNCUT\\\\n",dm1$myfill2,perl=T),]$myfill2 = gsub("UNCUT","1 UNCUT",dm1[grep("\\\\nUNCUT\\\\n",dm1$myfill2,perl=T),]$myfill2,perl=T)

dm2[grep("(^len$|total_read|_perkb_)",dm2$cat,invert=T,perl=T),]$mymean = dm2[grep("(^len$|total_read|_perkb_)",dm2$cat,invert=T,perl=T),]$mymean * 100
dm2[grep("(^len$|total_read|_perkb_)",dm2$cat,invert=T,perl=T),]$myse = dm2[grep("(^len$|total_read|_perkb_)",dm2$cat,invert=T,perl=T),]$myse * 100
# dm2[grep("\\\\nCUT\\\\n",dm2$myfill2,perl=T),]$myfill2 = gsub("CUT","2 CUT",dm2[grep("\\\\nCUT\\\\n",dm2$myfill2,perl=T),]$myfill2,perl=T)
# dm2[grep("\\\\nUNCUT\\\\n",dm2$myfill2,perl=T),]$myfill2 = gsub("UNCUT","1 UNCUT",dm2[grep("\\\\nUNCUT\\\\n",dm2$myfill2,perl=T),]$myfill2,perl=T)

dm1$myfill4 = dm1$myfill2
dm2$myfill4 = dm2$myfill2



dm1$myfill3 = dm1$myfill4
dm1$myfill3 = gsub("^([0-9]+)\\\\n(.+) (undigested|digested) (germline|enriched)\\\\n(UNCUT|CUT)\\\\nBAIT.+(PREY at.+)$","\\5MUHSPACEGROUP \\1 (\\2)MUHSPACE\\6",dm1$myfill3,perl=T)
dm1[grep("^UNCUT",dm1$myfill3),]$myfill3 = gsub("^UNCUT","A. NOT CUT",dm1[grep("^UNCUT",dm1$myfill3),]$myfill3,perl=T)
dm1[grep("^CUT",dm1$myfill3),]$myfill3 = gsub("^CUT","B. CUT",dm1[grep("^CUT",dm1$myfill3),]$myfill3,perl=T)
dm1$myfill3 = gsub("GROUP 2","GROUP 1",dm1$myfill3)
dm1$myfill3 = gsub("MUHSPACE","\\\n",dm1$myfill3,perl=T)
dm1$myfill5 = dm1$myfill3
dm1$myfill5 = gsub("\\\\n","\\\n",dm1$myfill5,perl=T)
dm1$myfill3 = gsub("^.+CUT\\\n","",dm1$myfill3,perl=T)
dm1$myfill2 = dm1$myfill3
dm1$myfill2 = gsub("\\\\n","\\\n",dm1$myfill2,perl=T)


dm2$myfill3 = dm2$myfill4
dm2$myfill3 = gsub("^([0-9]+)\\\\n(.+) (undigested|digested) (germline|enriched)\\\\n(UNCUT|CUT)\\\\nBAIT.+(PREY at.+)$","\\5MUHSPACEGROUP \\1 (\\2)MUHSPACE\\6",dm2$myfill3,perl=T)
dm2[grep("^UNCUT",dm2$myfill3),]$myfill3 = gsub("^UNCUT","A. NOT CUT",dm2[grep("^UNCUT",dm2$myfill3),]$myfill3,perl=T)
dm2[grep("^CUT",dm2$myfill3),]$myfill3 = gsub("^CUT","B. CUT",dm2[grep("^CUT",dm2$myfill3),]$myfill3,perl=T)
dm2$myfill3 = gsub("GROUP 2","GROUP 1",dm2$myfill3)
dm2$myfill3 = gsub("MUHSPACE","\\\n",dm2$myfill3,perl=T)
dm2$myfill5 = dm2$myfill3
dm2$myfill5 = gsub("\\\\n","\\\n",dm2$myfill5,perl=T)
dm2$myfill3 = gsub("^.+CUT\\\n","",dm2$myfill3,perl=T)
dm2$myfill2 = dm2$myfill3
dm2$myfill2 = gsub("\\\\n","\\\n",dm2$myfill2,perl=T)

# dm1$myfill2 = gsub("^(.+)([0-9]+)\\. (Ig[A-Z]+[0-9]+)  (germ|rich) (0[0-9]\\..+)$","\\1GROUP \\2 (\\3) \\4MUHSPACEPREY at \\5",dm1$myfill2,perl=T)
# dm1$myfill3 = dm1$myfill2
# dm1$myfill2 = gsub(" germ","",dm1$myfill2,perl=T)
# dm1$myfill2 = gsub(" rich","",dm1$myfill2,perl=T)
# dm1$myfill3 = gsub(" germ","",dm1$myfill3,perl=T)
# dm1$myfill3 = gsub(" rich","",dm1$myfill3,perl=T)
# dm1$myfill2 = gsub("GROUP 2","GROUP 1",dm1$myfill2,perl=T)
# dm1$myfill4 = gsub("(A|B)\\. (NOT CUT|CUT)","",dm1$myfill2,perl=T)
# dm1$myfill4 = gsub("GROUP 2","GROUP 1",dm1$myfill4,perl=T)
# dm1$myfill3 = gsub("MUHSPACE","\\\n",dm1$myfill3,perl=T)
# dm1$myfill4 = gsub("MUHSPACE","\\\n",dm1$myfill4,perl=T)

# get tech rep
dm5 = dm1
dm5$p = 1
dm5$pchar = ""
dm5$myse = 0

dmtemp1 = dm1


dmtemp1$p = 1
dmtemp1$pchar = ""

myfillz = unique(dmtemp1$myfill2)
cats = unique(dmtemp1$cat)
dms = dmtemp1
for (i in 1:length(myfillz)) {
   for (j in 1:length(cats)) {
      W = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "W" & dms$cat == cats[j],]$mymean
      S = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "S" & dms$cat == cats[j],]$mymean
      R = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "R" & dms$cat == cats[j],]$mymean
      D = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "D" & dms$cat == cats[j],]$mymean
      Wse = se(W)
      Sse = se(S)
      Rse = se(R)
      Dse = se(D)
      
      Wmean = mean(W)
      Smean = mean(S)
      Rmean = mean(R)
      Dmean = mean(D)
      
      print(paste(i,j,cats[j],length(W),length(S),length(R),length(D),Wmean,Smean,Rmean,Dmean))
      
      
      if (length(W) > 0) {
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "W" & dm5$cat == cats[j],]$mymean = Wmean
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "W" & dm5$cat == cats[j],]$myse = Wse
      }
      
      if (length(W) > 0 & length(S) > 0) {
         if (length(W) > 1 & length(S) > 1) {
            Sp = get_p(W,S)
            Spchar = get_pchar(Sp)
         } else {
            Sp = 1
            Spchar = ""
         }
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "S" & dm5$cat == cats[j],]$mymean = Smean
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "S" & dm5$cat == cats[j],]$myse = Sse
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "S" & dm5$cat == cats[j],]$p = Sp
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "S" & dm5$cat == cats[j],]$pchar = Spchar
      } else {
         Sp = 1
         Spchar = ""
      }
      if (length(W) > 0 & length(R) > 0) {
         if (length(W) > 1 & length(R) > 1) {
            Rp = get_p(W,R)
            Rpchar = get_pchar(Rp)
         } else {
            Rp = 1
            Rpchar = ""
         }
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "R" & dm5$cat == cats[j],]$mymean = Rmean
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "R" & dm5$cat == cats[j],]$myse = Rse
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "R" & dm5$cat == cats[j],]$p = Rp
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "R" & dm5$cat == cats[j],]$pchar = Rpchar
      } else {
         Rp = 1
         Rpchar = ""
      }
      if (length(W) > 0 & length(D) > 0) {
         if (length(W) > 1 & length(D) > 1) {
            Dp = get_p(W,D)
            Dpchar = get_pchar(Dp)
         } else {
            Dp = 1
            Dpchar = ""
         }
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "D" & dm5$cat == cats[j],]$mymean = Dmean
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "D" & dm5$cat == cats[j],]$myse = Dse
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "D" & dm5$cat == cats[j],]$p = Dp
         dm5[dm5$myfill2 == myfillz[i] & dm5$sampleType == "D" & dm5$cat == cats[j],]$pchar = Dpchar
      } else {
         Dp = 1
         Dpchar = ""
      }
      
      print(paste("    ",i,j,cats[j],length(W),length(S),length(R),length(D),Sp,Rp,Dp,Spchar,Rpchar,Dpchar))
      
      
   }
}



temp5 = dm5[grep("(IgG1.+IgG1|IgG3.+IgG3)",paste(dm5$sample,dm5$isotype),perl=T),]
temp5 = temp5[temp5$mypaste %in% dmtemp1$mypaste,]
temp5$p = as.integer(temp5$p * 1000)/1000

temp5[grep("_[A-Z]\\.\\.[A-Z]",temp5$cat),]$cat = gsub("\\.\\.","->",temp5[grep("_[A-Z]\\.\\.[A-Z]",temp5$cat),]$cat,perl=T)


temp5[grep("^mat_[A_Z]_any$",temp5$cat),]$pchar = ""
temp5[grep("^mat_any$",temp5$cat),]$pchar = ""


# get bio reps

dm = dm1
dm$p = 1
dm$pchar = ""
dm$myse = 0

dmtemp1 = dm1[dm1$mypaste %in% dm2$mypaste,]
dmtemp2 = dm2[dm2$mypaste %in% dm1$mypaste,]


dmtemp1$p = 1
dmtemp1$pchar = ""
dmtemp2$p = 1
dmtemp2$pchar = ""

myfillz = unique(dmtemp1$myfill2)
cats = unique(dmtemp1$cat)
dms = rbind(dmtemp1,dmtemp2)
for (i in 1:length(myfillz)) {
   for (j in 1:length(cats)) {
      W = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "W" & dms$cat == cats[j],]$mymean
      S = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "S" & dms$cat == cats[j],]$mymean
      R = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "R" & dms$cat == cats[j],]$mymean
      D = dms[dms$myfill2 == myfillz[i] & dms$sampleType == "D" & dms$cat == cats[j],]$mymean
      Wse = se(W)
      Sse = se(S)
      Rse = se(R)
      Dse = se(D)
      
      Wmean = mean(W)
      Smean = mean(S)
      Rmean = mean(R)
      Dmean = mean(D)

      print(paste(i,j,cats[j],length(W),length(S),length(R),length(D),Wmean,Smean,Rmean,Dmean))

      
      if (length(W) > 0) {
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "W" & dm$cat == cats[j],]$mymean = Wmean
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "W" & dm$cat == cats[j],]$myse = Wse
      }
      
      if (length(W) > 0 & length(S) > 0) {
         Sp = get_p(W,S)
         Spchar = get_pchar(Sp)
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "S" & dm$cat == cats[j],]$mymean = Smean
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "S" & dm$cat == cats[j],]$myse = Sse
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "S" & dm$cat == cats[j],]$p = Sp
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "S" & dm$cat == cats[j],]$pchar = Spchar
      } else {
         Sp = 1
         Spchar = ""
      }
      if (length(W) > 0 & length(R) > 0) {
         Rp = get_p(W,R)
         Rpchar = get_pchar(Rp)
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "R" & dm$cat == cats[j],]$mymean = Rmean
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "R" & dm$cat == cats[j],]$myse = Rse
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "R" & dm$cat == cats[j],]$p = Rp
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "R" & dm$cat == cats[j],]$pchar = Rpchar
      } else {
         Rp = 1
         Rpchar = ""
      }
      if (length(W) > 0 & length(D) > 0) {
         Dp = get_p(W,D)
         Dpchar = get_pchar(Dp)
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "D" & dm$cat == cats[j],]$mymean = Dmean
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "D" & dm$cat == cats[j],]$myse = Dse
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "D" & dm$cat == cats[j],]$p = Dp
         dm[dm$myfill2 == myfillz[i] & dm$sampleType == "D" & dm$cat == cats[j],]$pchar = Dpchar
      } else {
         Dp = 1
         Dpchar = ""
      }
      
      print(paste("    ",i,j,cats[j],length(W),length(S),length(R),length(D),Sp,Rp,Dp,Spchar,Rpchar,Dpchar))
      
      
   }
}

temp0 = dm[grep("(IgG1.+IgG1|IgG3.+IgG3)",paste(dm$sample,dm$isotype),perl=T),]
temp0 = temp0[temp0$mypaste %in% dmtemp1$mypaste,]
temp0$p = as.integer(temp0$p * 1000)/1000

dm1$myfill6 = dm1$myfill2
dm1$myfill2 = dm1$myfill5
dm2$myfill6 = dm2$myfill2
dm2$myfill2 = dm2$myfill5

temp1 = dm1[grep("(IgG1.+IgG1|IgG3.+IgG3)",paste(dm1$sample,dm1$isotype),perl=T),]
temp2 = dm2[grep("(IgG1.+IgG1|IgG3.+IgG3)",paste(dm2$sample,dm2$isotype),perl=T),]

temp0[grep("_[A-Z]\\.\\.[A-Z]",temp0$cat),]$cat = gsub("\\.\\.","->",temp0[grep("_[A-Z]\\.\\.[A-Z]",temp0$cat),]$cat,perl=T)
temp1[grep("_[A-Z]\\.\\.[A-Z]",temp1$cat),]$cat = gsub("\\.\\.","->",temp1[grep("_[A-Z]\\.\\.[A-Z]",temp1$cat),]$cat,perl=T)
temp2[grep("_[A-Z]\\.\\.[A-Z]",temp2$cat),]$cat = gsub("\\.\\.","->",temp2[grep("_[A-Z]\\.\\.[A-Z]",temp2$cat),]$cat,perl=T)

temp0[grep("^mat_[A_Z]_any$",temp0$cat),]$pchar = ""
temp0[grep("^mat_any$",temp0$cat),]$pchar = ""

#temp1$p = as.integer(temp1$p * 1000)/1000
#temp2$p = as.integer(temp2$p * 1000)/1000

p0.total = ggplot(temp0[temp0$cat == "total_read",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   #geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Combined; Total Read") +
   xlab("Sample") + ylab("Total read") +
   facet_grid(cat~myfill2,scales="free_y")

p0.len = ggplot(temp0[temp0$cat == "len",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
#   geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Combined; Average read length (bp)") +
   xlab("Sample") + ylab("average read length (bp)") +
   facet_grid(cat~myfill2,scales="free_y")

p0.perkb = ggplot(temp0[grep("_perkb_",temp0$cat),],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
 #  geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Combined; Mutation per kb") +
   xlab("Sample") + ylab("mutation per kb") +
   facet_grid(cat~myfill2,scales="free_y")

tempz0 = temp0
tempz0 = tempz0[grep("(total_read|^len$|_perkb_).+$",tempz0$cat,invert=T),]
tempz0 = tempz0[grep("(any|blunt)",tempz0$cat,perl=T),]
p0.any = ggplot(tempz0,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
  # geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Combined; % of read with at least 1 <mutation>, of total read)") +
   xlab("Sample") + ylab("% at least 1 <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")

tempz0 = temp0
tempz0 = tempz0[grep("^mat",tempz0$cat),]
tempz0 = tempz0[grep("(total_read|^len$|_perkb_).+$",tempz0$cat,invert=T),]
tempz0 = tempz0[grep("only",tempz0$cat),]
p0.only = ggplot(tempz0,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
#   geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Combined; % of read with at least 1 <mutation> ONLY, of total read)") +
   xlab("Sample") + ylab("% no <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")


#p5
p5.total = ggplot(temp5[temp5$cat == "total_read",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   #geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1 Combined; Total Read") +
   xlab("Sample") + ylab("Total read") +
   facet_grid(cat~myfill2,scales="free_y")

p5.len = ggplot(temp5[temp5$cat == "len",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   #geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1 Combined; Average read length (bp)") +
   xlab("Sample") + ylab("average read length (bp)") +
   facet_grid(cat~myfill2,scales="free_y")

p5.perkb = ggplot(temp5[grep("_perkb_",temp5$cat),],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   #geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1 Combined; Mutation per kb") +
   xlab("Sample") + ylab("mutation per kb") +
   facet_grid(cat~myfill2,scales="free_y")

tempz5 = temp5
tempz5 = tempz5[grep("(total_read|^len$|_perkb_).+$",tempz5$cat,invert=T),]
tempz5 = tempz5[grep("(any|blunt)",tempz5$cat,perl=T),]
p5.any = ggplot(tempz5,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   #geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1 Combined; % of read with at least 1 <mutation>, of total read)") +
   xlab("Sample") + ylab("% at least 1 <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")

tempz5 = temp5
tempz5 = tempz5[grep("^mat",tempz5$cat),]
tempz5 = tempz5[grep("(total_read|^len$|_perkb_).+$",tempz5$cat,invert=T),]
tempz5 = tempz5[grep("only",tempz5$cat),]
p5.only = ggplot(tempz5,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   geom_text(aes(group=treat,y=0,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   #   geom_text(aes(group=treat,y=-0.5,label=paste("p=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black",vjust=0) +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1 Combined; % of read with at least 1 <mutation> ONLY, of total read)") +
   xlab("Sample") + ylab("% no <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")


#p1 indiv
p1.total = ggplot(temp1[temp1$cat == "total_read",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
 #  geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1; Total Read") +
   xlab("Sample") + ylab("Total read") +
   facet_grid(cat~myfill2,scales="free_y")

p1.len = ggplot(temp1[temp1$cat == "len",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
#   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1; Average read length (bp)") +
   xlab("Sample") + ylab("average read length (bp)") +
   facet_grid(cat~myfill2,scales="free_y")

p1.perkb = ggplot(temp1[grep("_perkb_",temp1$cat),],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
#   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1; Mutation per kb") +
   xlab("Sample") + ylab("mutation per kb") +
   facet_grid(cat~myfill2,scales="free_y")


tempz1 = temp1
tempz1 = tempz1[grep("(total_read|^len$|_perkb_).+$",tempz1$cat,invert=T),]
tempz1 = tempz1[grep("(any|blunt)",tempz1$cat,perl=T),]
p1.any = ggplot(tempz1,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
 #  geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1; % of read with at least 1 <mutation>, of total read)") +
   xlab("Sample") + ylab("% at least 1 <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")

tempz1 = temp1
tempz1 = tempz1[grep("^mat",tempz1$cat),]
tempz1 = tempz1[grep("(total_read|^len$|_perkb_).+$",tempz1$cat,invert=T),]
tempz1 = tempz1[grep("only",tempz1$cat),]
p1.only = ggplot(tempz1,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
#   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep1; % of read with at least 1 <mutation> ONLY, of total read)") +
   xlab("Sample") + ylab("% no <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")


p2.total = ggplot(temp2[temp2$cat == "total_read",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
 #  geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep2; Total Read") +
   xlab("Sample") + ylab("Total read") +
   facet_grid(cat~myfill2,scales="free_y")

p2.len = ggplot(temp2[temp2$cat == "len",],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
 #  geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=2,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep2; Average read length (bp)") +
   xlab("Sample") + ylab("average read length (bp)") +
   facet_grid(cat~myfill2,scales="free_y")

p2.perkb = ggplot(temp2[grep("_perkb_",temp2$cat),],aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
#   geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep2; Mutation per kb") +
   xlab("Sample") + ylab("mutation per kb") +
   facet_grid(cat~myfill2,scales="free_y")


tempz2 = temp2
tempz2 = tempz2[grep("(total_read|^len$|_perkb_).+$",tempz2$cat,invert=T),]
tempz2 = tempz2[grep("(any|blune)",tempz2$cat,perl=T),]
p2.any = ggplot(tempz2,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
  # geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep2; % of read with at least 1 <mutation>, of total read)") +
   xlab("Sample") + ylab("% at least 1 <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")

tempz2 = temp2
tempz2 = tempz2[grep("^mat",tempz2$cat),]
tempz2 = tempz2[grep("(total_read|^len$|_perkb_).+$",tempz2$cat,invert=T),]
tempz2 = tempz2[grep("only",tempz2$cat),]
p2.only = ggplot(tempz2,aes(treat,mymean)) + 
   geom_bar(aes(fill=treat),stat="identity",position="dodge",color="black") +
 #  geom_errorbar(aes(group=treat,ymin=mymean-myse,ymax=mymean+myse),stat="identity",position=position_dodge(width=0.9),width=0.5,color="black") +
   geom_text(aes(color=treat,group=filedesc,label=as.integer(10*mymean)/10),stat="identity",position=position_dodge(width=0.9),color="black",size=3,vjust=0) +
   #geom_text(data=dmp,aes(group=treat,y=1,label=paste(pchar,"\np=",p,sep="")),stat="identity",position=position_dodge(width=0.9),size=3,color="black") +
   theme_bw() + theme(panel.grid=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank()) + ggtitle("Rep2; % of read with at least 1 <mutation> ONLY, of total read)") +
   xlab("Sample") + ylab("% no <mutation> of total read") +
   facet_grid(cat~myfill2,scales="free_y")


pdf("results_combined_IgG1and3.pdf",height=length(unique(temp1$cat)),width=10)
grid.arrange(p0.total,p0.len,p0.perkb,p0.any,p0.only,ncol=1,heights=c(2.3,2.3,16,50,11))
dev.off()


pdf("results_run1combiined_IgG1and3.pdf",height=length(unique(temp5$cat)),width=10)
grid.arrange(p5.total,p5.len,p5.perkb,p5.any,p5.only,ncol=1,heights=c(2.3,2.3,16,50,11))
dev.off()

pdf("results_run1_IgG1and3.pdf",height=length(unique(temp1$cat)),width=10)
grid.arrange(p1.total,p1.len,p1.perkb,p1.any,p1.only,ncol=1,heights=c(2.3,2.3,16,50,11))
dev.off()


pdf("results_run2_IgG1and3.pdf",height=length(unique(temp1$cat)),width=10)
grid.arrange(p2.total,p2.len,p2.perkb,p2.any,p2.only,ncol=1,heights=c(2.3,2.3,16,50,11))
dev.off()





