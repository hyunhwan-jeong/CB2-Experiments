pdf(file='sample1.pdf',width=4.5,height=4.5);
gstable=read.table('sample1.gene_summary.txt',header=T)
# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=3
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("PSMD6","PSMB2","NUP98","RPS13","PSMC4","PSMD11","PSMC1","RPL3","RPL11","COPB1")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='A1,A2,A3_vs_B1,B2,B3 neg.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5423.771421646316,5632.245093570229,5434.54446136728,673.157294769762,687.2494014427542,788.4373153160907),c(5046.296221815596,5174.490889581273,5020.756685733464,854.8189335382135,1004.9205749193144,999.6053881227292),c(4934.073865109166,4735.420530653091,5195.554863029849,1025.3790277152596,864.7715277973026,817.1286295561231),c(4229.206170093571,4162.760679540378,3901.850838971856,890.1420299654123,791.0635104220223,888.2830888714035),c(3291.546644224141,3092.1763575580458,3153.2801136008948,549.02184161132,551.2529186799131,394.792483942846),c(3668.0943865614195,3827.3856606994905,3693.4756106693853,948.6774469019134,827.3984485647661,709.2492880136012),c(5758.583576778724,5378.144800743707,5564.902424096787,2285.9089545030147,2092.8924370220434,2189.7211027992735),c(5172.43044092365,5111.900008627681,5229.131914035935,2180.9488965479095,1983.8876225938118,1876.4119512981197),c(4040.9322989249317,4133.800719696179,4100.350464037243,1406.8684691290077,1392.147201411984,1339.310548724713),c(3066.174473317839,3151.9646617525214,3017.9843492528453,758.9419575215308,953.0135204296804,892.8736991498087),c(3832.2543629005777,3690.0593995028034,3789.268962069099,1308.9730304593422,1383.8420726936426,1354.23003212953),c(1941.1685337732115,2040.2758806364868,2079.802041729873,456.1725595741116,376.8452155947429,436.10797644849265),c(2427.1562603365956,2518.5823141922933,2356.3189323682227,799.3112105811866,709.0503643284005,583.0075053574586),c(3797.010978149798,4087.0911070442444,3891.975235734772,1563.2993247351744,1785.6026744434098,1586.0558511889917),c(3900.886217415254,3722.7561283591576,4056.8978097940735,1810.5609997255667,1817.785048226983,1838.539416501277),c(2212.913579351592,2614.8041162552777,2314.8413987724703,1115.2006157729938,972.7382011357413,920.4173608202398),c(2135.0071499025003,2280.363289667429,2206.209763164547,1197.9575845452885,1147.1459042209117,1314.0621921934844),c(1980.1217484977574,2010.381728539249,1992.8967332435345,1163.643719444581,1173.0994314657287,1218.8070289165769),c(2246.302049115489,2236.4562537746106,2390.8835436980166,2613.9091356127187,2569.3991972368835,2416.9563115803303))
targetgene="PSMD6"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5802.174078970478,5714.454011837633,5639.957008698625,614.6218778332609,772.376970805754,740.2359073928362),c(5701.081212185347,5307.146189512767,5540.213416004078,724.6280924208232,603.1599731695471,574.9739373702496),c(4392.438688939288,4144.076834479604,4181.33041058133,212.9478098896848,166.10257436682883,164.1143174529854),c(4980.446739781245,5105.36066285641,5029.644728646839,599.48340793589,580.3208691941082,616.2894298758962),c(3963.953326969281,4145.011026732643,3649.0353961025075,236.16013039898695,288.60322296236507,183.62441113620744),c(4225.496340119805,4473.8466998022595,4314.651054281963,547.0033789583373,618.7320895164373,494.6382574981588),c(4155.937028111687,4007.6847655359566,4044.0595255858643,861.8835528236532,894.8776194012903,758.598348506457),c(3760.840135905577,3618.126596018825,3844.572340196769,565.1695428351825,675.8298494550348,755.1553907976531),c(3020.729056139202,3085.637011786775,2926.141239147965,301.76016662092775,379.95963886412096,254.77887045148782),c(3310.0957940929725,3269.6728856353957,3372.518505464158,613.6126465067696,579.2827281043155,454.4704175621134),c(5811.448653904894,5989.106534231006,6175.214704148573,2065.8965253278902,2093.930578111836,2257.43260440575),c(3411.1886608781037,3374.3024179757285,3395.232392909451,581.3172440590448,716.3173519569493,718.4305085704116),c(4005.688914174152,4289.810825953639,4029.2461207302385,1472.4685053509486,1215.6632161472285,1238.317122599799),c(2724.8701157313403,2710.0917260652236,2643.6989865673645,530.855677734475,508.6891339984133,577.2692425094522),c(1698.1746704915195,1667.5331716740518,1754.8946952298122,258.3632195817977,296.9083516807065,226.0875562114554),c(1387.4764101885926,1464.8134527646573,1563.307992430384,121.10775917896767,281.3362353338163,126.24178265614262),c(4433.246818650717,4581.2788089017085,4534.877006468935,2244.5304701168675,2188.4014172829698,2082.989413826353),c(3375.945276127324,3370.5656489635735,3257.9615079139844,2024.5180409417428,1822.9757536759464,1964.7811991574195),c(715.9971849368941,747.3538024309476,680.4290630350819,229.09551111354716,180.63654962392636,221.49694593305023),c(6530.228211322113,6338.494436867474,6533.69910165472,4850.365755117655,4583.392911434683,5006.060508600855),c(3769.1872533465507,3838.5959677359547,3807.04504789585,3324.4079894626625,3820.359210437063,3224.903720579643),c(2467.036932554583,2660.5795366541734,2522.2290667512325,2846.0323407057404,3019.9524302069067,2917.3328319264956),c(3417.680863332195,3467.7216432795967,3448.560650389704,4063.1653204543654,4257.416609239782,3669.0452650153447),c(3058.7548133703062,2927.7585210232373,3328.07829089728,4149.959214532625,4069.5130719873064,3627.7297725096983))
targetgene="PSMB2"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(4819.996593415853,5083.874241036521,4684.98617567261,1257.502232808281,998.6917283805583,1271.5990471182365),c(3497.44220776817,3634.0078643204824,3323.140489278738,466.26487283902554,703.8596588794371,866.4776900489788),c(4315.459716983637,4108.577528864134,4487.4741109309325,1229.2437556665218,1242.654884481838,1203.88754551176),c(4931.291492628841,5051.177512180167,4897.311645269915,1532.013153613941,1828.1664591249098,1843.130026779682),c(6689.750900194063,6884.996904895104,6147.563015084738,3274.955654464584,3098.8511530311503,2935.6952730401163),c(3931.4923146988262,4110.445913370212,4058.87293044149,1400.8130811700594,1684.9029887335198,1294.5520985102623),c(3099.5629430817357,2994.0861709889837,3011.0714269868863,870.9666347620758,1006.9968570988998,1133.8807387660809),c(4735.59796151267,4823.234602438728,4888.423602356539,2726.9430441797554,2974.2742222560287,3011.4403426338017),c(3511.354070169794,3580.7589058972776,3500.901347546249,2207.1889110366856,2332.7030287641524,2040.5262687511051),c(2649.746058762573,2839.010256984562,3135.5040277741437,1646.0562935074688,1749.267736300666,1677.8680567570955),c(4330.299036878702,4418.729356872977,4352.178346582882,3095.3124783491153,2938.9774252030775,2646.4868255005895),c(7990.046305999148,7923.818690274122,7942.947683586595,5989.787922726443,6005.646204450654,5735.96754286728),c(6745.398349800557,6536.543194511675,6565.301032013388,5442.784543768105,5112.84486722895,4752.429290718968),c(2090.4891902173044,2139.3002594585873,2232.8738919046737,1393.7484618846197,1231.2353324941187,1665.2438784914812),c(4359.050219175391,4414.058395607784,3976.9054235736935,3340.555690686525,3591.968170682673,3568.0518388904306),c(4074.320768688828,4047.85503241662,3857.410624404978,3365.7864738488097,3813.092222808514,3877.918032682781))
targetgene="NUP98"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5344.937534703783,5181.030235352544,5067.172020947758,188.72625805389129,198.2849481504019,224.9399036418541),c(4575.1478151472775,4706.460570808892,4436.120974098096,649.9449742604598,645.723757851047,757.4506959368557),c(4966.534877379621,4740.091491918285,5083.9605464508,951.7051408813876,761.9955599078272,1028.2967023627616),c(4882.1362454764385,4746.630837689556,4692.886658262278,1179.7914206684434,985.1958942132535,1180.9344941197342),c(3169.1222550898533,3289.290922949208,3029.8350731373457,446.0802463091976,250.19200264003592,335.11455032357856),c(4406.350551340911,4704.592186302815,5008.905961848963,1163.643719444581,1045.4080774212289,1273.894352257439),c(4474.054948362146,4570.068501865245,4458.834861543389,1127.3113916908908,1229.1590503145333,996.1624304139253),c(2422.5189728693877,2475.609470552514,2714.803329874369,402.6832992700675,451.59137405981585,591.0410733446677),c(1537.7245241261273,1733.8608216397984,1668.976947067182,505.62489457219004,624.9609360551934,589.8934207750664),c(2784.227395311601,2717.5652640895332,2981.4446172756348,1618.8070476922012,1391.1090603221915,1685.9016247443044))
targetgene="RPS13"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(10275.301569839183,10376.073354500668,9832.150582840748,3859.300592503103,3662.5617647885756,3514.11216811917),c(4651.199329609486,4437.413201933751,4368.966872085925,970.8805360847241,1167.9087260167653,1044.3638383371797),c(3328.644943961804,3275.2780391536276,3492.0133046328733,471.3110294714825,424.5997057252062,488.8999946501523),c(2540.3060745364673,2581.173195145885,2348.4184497785554,237.16936172547835,260.5734135379627,261.66478586909557),c(2492.0782848775057,2219.6407932199145,2544.9429541965255,369.3786654958514,303.1371982194626,378.72534796842785),c(3915.7255373103194,3846.069505760264,3909.751321561523,1503.754676472182,1392.147201411984,1346.1964641423208),c(2971.5738089867987,3091.242165305007,3117.7279419473925,1472.4685053509486,1653.7587560397394,1620.4854282770307),c(3998.2692542266195,3718.085167093964,3632.2468705994647,2529.1337041874417,2376.3049545354447,2472.0436349211927),c(1713.941447880026,1860.9109680530594,2022.5235429547863,1042.5359602656133,1075.5141690252167,869.9206477577827),c(3192.308692425893,3298.632845479595,3143.404510363811,2210.21660501616,2121.960387536238,2045.1168790295103))
targetgene="PSMC4"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(8725.520098298315,8943.022438339327,9267.266077679547,1023.3605650622768,1241.6167433920455,914.6790979722333),c(3852.6584277562924,3885.3055803878888,3710.2641361724277,534.8926030404406,351.92982943971856,286.9131424003241),c(5254.97415783995,5362.263532442049,5671.558939057294,1364.480753416369,1225.0064859553627,1145.3572644620938),c(2382.6383006513997,2522.319083204448,2405.6969485536424,590.4003259974673,562.6724706676326,526.7725294469951),c(2775.8802778706267,2901.601137938154,2821.459844834875,934.5482083310338,930.1744164542414,727.6117291272219),c(2772.1704478968604,3072.558320244233,2950.8302472406745,986.0190059820951,1074.4760279354239,965.1758110346904),c(2050.608517999317,2176.667949580135,2257.5628999973837,872.9850974150586,875.1529386952294,704.6586777351961),c(1610.06620861457,1732.9266293867597,1671.9396280383073,545.9941476318459,516.9942627167547,546.2826231302171),c(2640.471483828157,2647.500845111632,2554.8185574336094,1216.1237484221338,1276.9135404449967,1253.2366060046159),c(2249.0844215958136,2110.340299614388,2087.7025243195403,1389.711536578654,1176.2138547351067,1279.6326151054457))
targetgene="PSMD11"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(6292.79909300107,6266.561633383495,6171.26446285374,869.9574034355844,960.2805080582292,765.4842639240647),c(6955.0037433183525,7082.111470286267,7154.874545267298,1350.3515148454894,1359.964827628411,1343.9011590031182),c(7253.645056206539,7230.648038519418,7001.802695092498,1393.7484618846197,1517.7622732768984,1357.6729898383337),c(6400.384162240292,6069.447067992333,6672.9451072976035,1157.5883314856326,1230.197191404326,1062.7262794508006),c(4750.437281407735,4835.379101728231,4623.75743560269,963.8159167992843,1237.4641790328747,1078.7934154252187),c(4374.816996563898,4303.823709749219,4237.621349032709,1434.1177149442756,1532.2962485339958,1339.310548724713),c(3852.6584277562924,3779.7418557945175,4079.6116972393665,2355.545916030921,2577.7043259552247,2274.6473929497697),c(4542.686802876822,4705.526378555854,4195.156255113248,2807.681550299067,3111.3088461086622,3273.1051285028975),c(5548.978183260929,5235.213386028788,5249.87068083381,4969.45505164364,4446.358287582049,4777.677647250197),c(3606.8821919942757,3567.680214354736,3668.7866025766753,4121.700737390866,3857.7322896695996,4183.193616196726))
targetgene="PSMC1"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5249.4094128793,5478.1033718188455,5018.781565086047,1268.6037773996864,1193.8622532615823,1254.384258574217),c(2764.750787949328,2642.8298838464384,2591.35828941082,331.0278750891783,333.2432898234503,336.26220289317985),c(5075.974861605727,5176.359274087351,5138.2763642547625,1960.9364673727848,1960.0103775285802,2038.2309636119026),c(1248.3577861723568,1367.657458448634,1199.8857933056959,293.6863160089966,400.72246065997456,262.8124384386969),c(4736.525419006111,5005.402091781271,4892.373843651373,2825.847714175912,3005.418454949809,2924.2187473441036),c(2905.724326952447,2760.5381077293127,2700.9774853424515,1698.5363224850216,1799.0985086107148,1768.5326097555978),c(4257.95735239026,4471.978315296183,3992.706388753028,3161.9217458975477,3013.7235836681507,2831.258889206398),c(3683.8611639499263,3909.5945789668945,3693.4756106693853,2771.349222545377,2621.3062517265175,2592.5471547293287),c(2401.1874505202313,2346.6909396331753,2629.873142035447,1790.3763731957388,1545.7920827013008,1544.740358683345),c(7402.038255157191,7268.015728640965,7688.15712006983,5725.369315185697,6038.866719324021,5943.692657965114),c(5681.604604823074,5732.203664645368,5826.605909879511,4783.756487569223,5158.523075179828,5923.034911712291))
targetgene="RPL3"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(6121.219456714379,6392.677587543718,5979.677760054312,1017.3051771033284,891.7631961319123,868.7729951881814),c(3581.8408396713535,3497.6157953768347,3411.033358088785,121.10775917896767,91.35641590175585,175.59084314899835),c(6660.9997178973745,6678.540416973556,6387.540173745878,1886.2533492124214,1955.8578131694094,2053.150447016719),c(3522.483560091093,3448.1036059657845,3367.580703845616,297.7232413149622,440.1718220720964,420.0408404740745),c(3540.1052524664824,3647.0865558630244,3776.4306778608902,931.5205143515597,860.6189634381318,829.7528078217373),c(2762.8958729624446,2496.161700119365,2582.4702464974444,418.8310004939299,369.57822796619416,285.7654898307228),c(2276.9081463990606,2340.1515938619045,2447.174482149395,385.52636671971374,240.8487328319018,242.15469218587356),c(2256.5040815433463,2040.2758806364868,2199.2968408985885,375.43405345479977,553.3292008594985,485.4570369413484),c(5753.946289311516,5516.405254193432,5758.464247543632,3201.2817676307122,3579.510477605161,3849.2267184427483),c(2884.3928046032906,2967.9287879039007,2970.5814537148426,2084.0626892047353,2374.2286723558595,2117.418990914392))
targetgene="RPL11"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5770.640524193464,5717.256588596749,5900.672934157641,561.1326175292169,827.3984485647661,719.5781611400129),c(6983.754925615041,7000.836744271902,6491.234007735259,2392.887475111103,2218.5075088869576,2242.513121000933),c(3721.8869211810306,4102.038183092864,3989.7437077819027,791.2373599692554,776.5295351649247,867.6253426185801),c(4846.8928607256585,5105.36066285641,4525.988963555559,1405.8592378025164,1552.0209292400568,1506.8678238865023),c(4051.134331352789,3828.319852952529,4044.0595255858643,1145.4775555677359,940.5558273521682,1175.1962312717276),c(2521.7569246676358,2298.112942475164,2603.209013295321,354.24019559848045,371.6545101457795,228.382861350658),c(3190.4537774390096,3202.4110434166105,3175.994001046188,916.3820444541886,813.9026143974612,967.4711161738929),c(2799.0667152066662,2530.7268134817964,2931.0790407665068,922.4374324131371,980.0051887642901,1078.7934154252187),c(1844.712954455288,1785.2413955569261,1807.2353923863568,479.3848800834137,622.8846538756081,548.5779282694197),c(6751.890552254648,6536.543194511675,6316.435830438873,7781.173527248673,7526.522900996931,7579.097569646962))
targetgene="COPB1"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



# 
#
# parameters
# Do not modify the variables beginning with "__"

# gstablename='__GENE_SUMMARY_FILE__'
startindex=9
# outputfile='__OUTPUT_FILE__'
targetgenelist=c("OR52E8","FAM71B","TGM6","TAAR1","TAS2R9","NLRP5","DEFB129","PLA2G2E","KRT77","VN1R2")
# samplelabel=sub('.\\w+.\\w+$','',colnames(gstable)[startindex]);
samplelabel='A1,A2,A3_vs_B1,B2,B3 pos.'


# You need to write some codes in front of this code:
# gstable=read.table(gstablename,header=T)
# pdf(file=outputfile,width=6,height=6)


# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")

######
# function definition

plotrankedvalues<-function(val, tglist, ...){
  
  plot(val,log='y',ylim=c(max(val),min(val)),type='l',lwd=2, ...)
  if(length(tglist)>0){
    for(i in 1:length(tglist)){
      targetgene=tglist[i];
      tx=which(names(val)==targetgene);ty=val[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      # text(tx+50,ty,targetgene,col=colors[i])
    }
    legend('topright',tglist,pch=20,pt.cex = 2,cex=1,col=colors)
  }
}



plotrandvalues<-function(val,targetgenelist, ...){
  # choose the one with the best distance distribution
  
  mindiffvalue=0;
  randval=val;
  for(i in 1:20){
    randval0=sample(val)
    vindex=sort(which(names(randval0) %in% targetgenelist))
    if(max(vindex)>0.9*length(val)){
      # print('pass...')
      next;
    }
    mindiffind=min(diff(vindex));
    if (mindiffind > mindiffvalue){
      mindiffvalue=mindiffind;
      randval=randval0;
      # print(paste('Diff: ',mindiffvalue))
    }
  }
  plot(randval,log='y',ylim=c(max(randval),min(randval)),pch=20,col='grey', ...)
  
  if(length(targetgenelist)>0){
    for(i in 1:length(targetgenelist)){
      targetgene=targetgenelist[i];
      tx=which(names(randval)==targetgene);ty=randval[targetgene];
      points(tx,ty,col=colors[(i %% length(colors)) ],cex=2,pch=20)
      text(tx+50,ty,targetgene,col=colors[i])
    }
  }
  
}




# set.seed(1235)



pvec=gstable[,startindex]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='RRA score',main=paste('Distribution of RRA scores in \\n',samplelabel))


pvec=gstable[,startindex+1]
names(pvec)=gstable[,'id']
pvec=sort(pvec);

plotrankedvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))

# plotrandvalues(pvec,targetgenelist,xlab='Genes',ylab='p value',main=paste('Distribution of p values in \\n',samplelabel))



# you need to write after this code:
# dev.off()






# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(4340.50106930656,4203.86513867408,4491.424352225766,5401.406059381958,5274.794877236608,5554.638436870275),c(4202.309902783765,4281.403095676291,4210.957220292583,5582.058466823918,5804.246833030875,5714.162144044855),c(4320.097004450845,4644.803882108339,4329.464459137589,5555.818452335142,6145.7952515726665,5646.450642438379),c(1804.8322822373002,1868.384506077369,1869.4516927799855,2845.023109379249,2435.4789966536277,2245.956078709737),c(7547.649081627517,7749.124738955888,8007.13910462764,9602.836071565645,10719.844893199215,9821.610690647894),c(6658.21734541705,6692.553300769136,6712.447520245939,8002.195187750289,8063.2418444197465,9304.019381757711),c(4726.323386578254,4698.987032784583,4749.177596713656,6118.969532517342,6430.245910175861,6209.948054112615),c(6353.083830074772,6408.558855845376,6606.778565609141,8686.454027111457,8604.113352201734,8180.467516118041),c(8717.17298085734,8665.567339186837,8676.70500410193,11286.233924153295,10959.655484941324,10636.444015064815),c(5634.304272657553,5605.153518232107,5248.883120510102,7755.942744086387,7262.83506418959,6992.647106580699))
targetgene="OR52E8"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(5357.921939611965,5649.994746377964,5660.6957754965015,5277.270606223517,5758.568625079997,5803.679044473756),c(3517.846272623885,3582.627290403355,3443.622848771162,3593.8727536358656,4031.101851664977,3845.7837607339443),c(6596.077693356464,6251.614557334877,6302.609985906956,7489.505673892659,7124.762299247164,7481.547101230852),c(4725.395929084812,4451.426085729331,4633.633038839774,5878.772476812389,6161.367367919557,6040.095473811623),c(4669.748479478318,4769.051451762484,4693.874218585986,7093.886993908031,6814.358113399153,6367.176456147993),c(8398.12760311344,8685.185376500649,8320.195727243201,10216.448718072414,10772.790088778642,10059.174772555363),c(7091.339994854264,7217.569346976877,7118.334813290087,9100.23887097293,9207.27332537128,10056.879467416162),c(4110.49161093305,3864.7533508210377,3864.323546670937,5003.768916744348,5309.053533199766,5233.295717381912),c(9619.58912197599,9980.910031465304,9864.740073523126,11804.978825969874,12454.578654242783,12284.473105012277),c(5837.417463721258,6004.987802532663,5570.827786039038,7828.607399593769,7862.88061408976,8164.400380143623))
targetgene="FAM71B"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2071.0125828550317,1982.3559609480885,2258.550460321092,2036.6288168596395,2308.8257836989205,2419.251616719533),c(2381.7108431579586,2390.5979755259937,2549.8807558150675,2684.555328467117,2955.5876826397603,2468.6006772123887),c(2774.0253628837436,2837.141872478485,2797.758397065874,3597.9096789418313,3143.4912198922357,3119.319684176324),c(4692.934916814357,4582.213001154747,4586.230143301771,5622.427719883574,5212.506411849047,5669.403693830404),c(4639.1423821947465,4613.975537758062,4497.349714168016,5901.984797321691,5647.48752847218,5773.840077664123),c(3629.1411718368736,3534.0492932453435,4023.320758787988,4573.836371659012,4278.179431035635,4358.784459345724),c(5126.057566251571,4935.33767280337,4755.102958655906,6365.221976181242,6283.868016515093,6234.048758074243),c(4337.718696826235,4160.8922950343,4295.887408131504,5462.9691702979335,5707.699711680156,5660.222473273594),c(4577.002730134161,4424.33451039121,4833.120224228869,5496.27380407215,6052.362553491325,5197.718487724272),c(6941.091880916729,6913.9568647393035,6835.892560709488,8858.023352614993,8393.370710973819,9244.341448138443),c(3838.7465653546687,4135.669104202256,3802.107246277308,4864.494993688535,5197.972436591949,5183.946656889056),c(3333.282231429012,3485.4712960873317,3532.5032779049175,4590.993304209366,4488.922072263549,4387.4757735857565),c(5101.943671422091,5029.6910903602775,5068.159581271467,7010.120793809245,6398.0635363922875,6667.861429383533),c(4264.449554844351,4390.703589281817,4422.295129566178,5337.824485813,5052.632684020974,5722.195712032064))
targetgene="TGM6"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3653.255066666354,3726.492897371312,3549.29180340796,4233.725414631412,4259.492891419367,3961.6966702636755),c(1568.3306214096992,1524.6017569591331,1521.8304588346316,2016.4441903298116,1619.500100076581,1999.2107762454584),c(3225.6971621897897,3405.130762326005,3636.1971118942984,3945.0852552548718,4102.733586860672,3882.508642961186),c(3268.360206888102,3535.917677751421,3207.5959314048564,4440.617836562148,3658.409200429405,4042.032350135766),c(4485.184438283445,4812.958487655303,4461.797542514514,6350.083506283871,6227.808397666288,6083.706271456473),c(5099.161298941766,5221.200502233208,5132.351002312512,6212.828045881041,5983.845241565008,6075.672703469263),c(3687.5709939236926,3628.4027108022506,3709.2765758487194,4553.651745129185,4774.410871956536,4325.502534827287),c(4849.675233205983,4682.171572229887,4982.2418331088365,6400.545072608441,6378.338855686226,6492.270586234534),c(8873.913297248966,9048.586162932697,8355.747898896703,10636.288949892836,11063.469593920592,11303.230158003169),c(2437.358292764453,2414.8869741049994,2548.893195491359,3334.5003027275766,3351.1194378507716,3332.783062122165))
targetgene="TAAR1"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(7056.096610103484,7453.919986995664,7042.292668364542,7901.272055101149,8345.616220843356,7697.305784315896),c(3664.384556587653,3650.823324875179,3566.0803289110027,4452.728612480045,3896.1435099919286,4341.569670801705),c(5978.391002724376,6123.630218668577,6339.149717884166,6821.394535755354,6645.141115762945,7835.024092668051),c(4713.338981670072,4439.281586439829,4395.631000826052,6252.1880676142055,5994.226652462935,5921.88725914269),c(3303.6035916388814,3281.8173849248988,3340.9165751054893,4416.396284726355,3615.845415747905,3797.58235281069),c(4924.79929017475,5053.980088939283,5079.0227448322585,6398.526609955458,6328.508083376179,6475.055797690515),c(3519.701187610768,3564.87763759562,3581.8812940903367,4645.491795839902,3983.3473615345138,4231.39502411998),c(3402.84154343713,3511.628679172415,3574.968371824378,4370.980875034242,4764.029461058609,4750.133985579766),c(2447.56032519231,2769.8800302596997,2463.963007652437,3329.4541460951195,3560.823937988893,3378.689164906217),c(5179.850100871183,5010.073053046465,5072.1098225663,6318.797335162638,5974.501971756874,6115.840543405309))
targetgene="TAS2R9"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2877.9006021491996,2662.4479211602506,2668.3879946600746,2879.3369744799566,2967.00723462748,3122.762641885128),c(4619.665774832473,4518.687927948116,4390.69319920751,5027.990468580141,4766.105743238195,5156.402995218625),c(2822.2531525427053,3064.150589966885,2832.3230083956673,3965.2698817847,3364.615272018076,3223.756068010042),c(6446.757036912371,6470.215544545928,6180.152505767115,8787.377159760596,8565.702131879403,8614.280187427332),c(4413.7702112884435,4111.3801056232505,4066.7734130311574,5034.045856539089,5272.718595057023,5536.275995756654),c(6243.643845848666,5970.422689170233,6033.993577858274,7973.936710608529,8452.544753092001,8310.152256482988),c(6159.245213945484,6129.2353721868085,6085.34671469111,9248.595875967165,8776.444773107318,8792.166335715532),c(5194.689420766248,4904.509328453093,5139.26392457847,7401.702548487908,7153.830249761359,6901.982553582197),c(7210.982011508227,6767.28868101223,6994.889772826539,8431.118501509132,8136.949861795027,8765.770326614702),c(3097.7080280948526,3116.4653561370515,3204.633250433731,3957.1960311727685,4241.844492892891,4177.455353348719))
targetgene="NLRP5"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3255.3758019799197,3310.777344769098,3170.0686391039376,3723.054363426764,3260.801163038808,3809.058878506703),c(4415.625126275327,4426.2028948972875,4525.988963555559,5144.052071126652,4755.724332340268,5314.779049823604),c(3019.8015986457604,3002.4939012663317,2996.2580221312605,3567.632739147089,3426.903737405637,3573.790101738437),c(1677.7706056358047,1800.188471605545,1903.028743786071,2349.4905280719727,2113.655258817897,2228.7412901657176),c(4095.6522910379845,4106.709144358057,4043.071965262156,5810.144746610974,4853.30959478078,4678.979526264486),c(6427.280429550098,6202.102367923826,6489.258887087843,7087.831605949083,7523.408477727553,8186.205778966048),c(8555.795376998507,8303.100745007827,7949.860605852554,10472.79347500123,10590.07725697513,10986.478048793211),c(6525.590923854905,6664.527533177975,6517.898136475385,9628.06685472793,8544.93931008355,9570.274777905212),c(4331.226494372144,4470.109930790105,4556.6033335905195,5474.070714889339,5574.817652186693,5221.819191685899),c(6481.072964169709,6156.326947524931,6050.782103361316,8866.097203226926,8830.428109776538,8778.394504880316))
targetgene="DEFB129"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2516.1921797069863,2601.725424712736,2681.2262788682838,2994.389345699976,2741.7306181424683,3325.897146704557),c(6315.98553033711,5881.674425131558,6182.127626414533,7026.268495033108,7312.665836499639,7243.983019323383),c(4539.904430396497,4646.672266614417,4486.486550607224,5362.046037648794,5491.766365003278,5318.222007532408),c(4676.240681932409,4325.310131569109,4582.279902006938,5718.304695900257,5616.3432957784,5135.745248965802),c(8834.96008252442,8635.6731870896,8919.644843734195,11413.39707129121,11441.352950605127,10896.96114836431),c(7561.560944029141,7074.637932261958,7542.985752484697,8173.764513253826,8735.957270605404,9335.006001136946),c(2919.6361893540707,2917.4824062398116,2798.745957389582,3346.6110786454733,3386.416234903723,3630.0250776489006),c(4873.789128035464,4613.975537758062,4812.381457430993,6083.646436090143,5959.967996499777,5558.081394579079),c(4195.817700329674,4343.059784376844,4187.255772523581,4656.593340431307,6076.239798556557,5170.174826053841),c(4592.7695075226675,4471.978315296183,4641.533521429441,5133.959757861738,5713.928558218911,5828.927401004985))
targetgene="PLA2G2E"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(2801.849087686991,2592.3835021823493,2592.345849734528,2958.057017946285,3070.821343606748,3084.890107088285),c(4115.128898400258,3968.4486909083316,3841.609659225644,4353.823942483888,4580.278488165305,4715.704408491727),c(2170.25053465328,2126.221567916046,2213.122685430506,2465.5521306184833,2514.3777194778713,2878.3126445600515),c(963.6283356857939,1143.4513177193498,971.7593585290575,1558.2531681027174,1189.7096889024115,1262.4178265614262),c(5420.989049165992,5591.140634436527,5567.865105067913,6385.40660271107,6502.915786461349,6325.860963642346),c(6720.356997477635,7024.191550597869,7020.566341242957,8532.041634158273,9408.67269679106,8520.172676720025),c(7239.733193804916,7043.809587911681,7080.807520989169,8531.032402831781,8626.952456177172,8578.70295776969),c(7328.769113175306,7378.25041449953,7147.96162300134,8604.706289665653,8071.546973138088,8234.407186889302),c(4568.655612693186,4853.1287545359655,4782.754647719741,6048.323339662944,5316.320520828315,5668.256041260804),c(4427.682073690067,4793.34045034149,4508.212877728808,6134.108002414712,5767.9118948881305,5530.537732908648))
targetgene="KRT77"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}




# parameters
# Do not modify the variables beginning with "__"
targetmat=list(c(3339.774433883103,3481.734527075177,3613.4832244490053,3796.7282502606363,3485.0396384340274,4085.6431477806154),c(5524.864288431448,4928.798327032099,5345.664032233524,5488.199953460218,6325.3936601068,6068.786788051656),c(3252.593429499595,3044.532552653073,3233.272499821275,3632.223544042539,3542.137398372625,4070.7236643757988),c(3211.785299788166,3094.9789343171615,3246.1107840294835,3942.057561275398,3868.1137005675264,3732.1661563434163),c(3569.783892256613,3514.431255931531,3729.027782322887,4830.1811285878275,4819.0509388176215,4672.093610846878),c(4597.406794989875,4683.105764482925,4766.953682540407,6363.2035135282595,5886.259979124497,5750.887026272097),c(5947.784905440805,6016.198109569128,5729.824998156089,7595.474963174255,7270.102051818139,7420.721515041983),c(4822.778965896177,4974.573747430994,5059.271538358091,6127.043383129273,5825.009654826728,6207.652748973413),c(5506.315138562616,5380.013185249784,5340.726230614982,6711.388321167791,6617.111306338544,6649.498988269912),c(4821.851508402736,4907.311905212209,4594.130625891438,6022.083325174167,6571.433098387665,5854.1757575362135))
targetgene="VN1R2"
collabel=c("B1","B2","B3","A1","A2","A3")

# set up color using RColorBrewer
#library(RColorBrewer)
#colors <- brewer.pal(length(targetgenelist), "Set1")

colors=c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00",  "#A65628", "#F781BF",
          "#999999", "#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3", 
          "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5",
          "#D9D9D9", "#BC80BD", "#CCEBC5", "#FFED6F")


## code

targetmatvec=unlist(targetmat)+1
yrange=range(targetmatvec[targetmatvec>0]);
# yrange[1]=1; # set the minimum value to 1
for(i in 1:length(targetmat)){
  vali=targetmat[[i]]+1;
  if(i==1){
    plot(1:length(vali),vali,type='b',las=1,pch=20,main=paste('sgRNAs in',targetgene),ylab='Read counts',xlab='Samples',xlim=c(0.7,length(vali)+0.3),ylim = yrange,col=colors[(i %% length(colors))],xaxt='n',log='y')
    axis(1,at=1:length(vali),labels=(collabel),las=2)
    # lines(0:100,rep(1,101),col='black');
  }else{
    lines(1:length(vali),vali,type='b',pch=20,col=colors[(i %% length(colors))])
  }
}



dev.off()
Sweave("sample1_summary.Rnw");
library(tools);

texi2dvi("sample1_summary.tex",pdf=TRUE);

