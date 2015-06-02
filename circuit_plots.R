# generate circuit plots - using some code from HG #

setwd('~/stochastic_flies/gillespie')

load_gillespie = function(filename,ntrials) {
  con = file(filename,'r') # open connection 
  # using scan to retrieve header info line by line
  times = head(scan(con,nlines=1,sep='\t'),-1)
  species = head(scan(con,nlines=1,sep='\t',what=character()),-1)
  reactions = head(scan(con,nlines=1,sep='\t',what=character()),-1)
  parameters = head(scan(con,nlines=1,sep='\t',what=character()),-1) # could be parsed out
  
  concentrations.mat = matrix(nrow=ntrials*length(times),ncol=length(species)+2)
  colnames(concentrations.mat) = c(species,'times','trial')
  rxn_counts.mat = matrix(nrow=ntrials,ncol=length(reactions))
  colnames(rxn_counts.mat) = reactions
  for (i in 1:ntrials) {
    # concentration data in times x species format
    rdata = matrix(scan(con,nlines=length(times),sep='\t'),
                   nrow=length(times),ncol=length(species)+1,byrow=T)
    concentrations.mat[(1+(i-1)*length(times)):(i*length(times)),] = cbind(rdata[,1:length(species)],times,i)
    
    # reaction counts in one line 
    rxndata = head(scan(con,nlines=1,sep='\t'),-1)
    rxn_counts.mat[i,] = rxndata
  }
  
  close(con)
  
  ## return extracted data ##
  results = list()
  results$concentrations.mat = concentrations.mat
  results$rxn_counts.mat = rxn_counts.mat
  results$parameters = parameters
  return(results)
}

pl300 = load_gillespie('results/yan_network_pl300.txt',100)
ERKpl300 = load_gillespie('results/yan_network_ERK_pl300.txt', 100)

ERKsub300 = load_gillespie('results/yan_subnetwork.txt', 100)

## plotting library things ##
require(ggplot2)
require(reshape2)
require(plyr)
require(ggthemes)

convert_concs = function(concentrations.mat) {
  # get total Yan column
  bonus_concs.df = data.frame(concentrations.mat)
  bonus_concs.df$YTotal = with(bonus_concs.df,Y + YP + M_Y + M_YP + 2*Y_Y)
  return(bonus_concs.df)
}

plot_circuit1_course = function(concentrations.mat,outfile,plot_legend=T,trial=1) {
  bonus_concs.df = convert_concs(concentrations.mat)
  species_toplot = c('YTotal','M','P2P')
  species_names = c('Total Yan','Mae','pPntP2')
  sub_conc.df = bonus_concs.df[,c(species_toplot,'times','trial')]
  sub_conc.mlt = melt(sub_conc.df,id.vars=c('times','trial'))
  sub_conc.summary = ddply(sub_conc.mlt[,c('times','variable','value')],c('times','variable'),
                           summarise,mean=mean(value))
  
  # overlay single stochastic run
  trial_sub_conc.mlt = sub_conc.mlt[sub_conc.mlt$trial %in% trial,]
  
  # plot of these values
  g = ggplot(sub_conc.summary,aes(times/3600,mean,col=variable))
  g = g + geom_line(size=1.5)
  # plot stochastic trial overlay
  g = g + geom_line(data=trial_sub_conc.mlt,aes(times/3600,value,col=variable),size=0.7,alpha=0.5)
  ## axis ##
  g = g + scale_x_continuous(breaks=seq(0,60,10),limits=c(0,60))
  g = g + scale_y_continuous(breaks=seq(0,800,100),limits=c(0,800))
  ## labels ##
  if (plot_legend) {
    g = g + scale_colour_discrete(breaks=species_toplot,labels=species_names)
  }
  else {
    g = g + scale_colour_discrete(guide=F,breaks=species_toplot,labels=species_names)
  }
  ## theme things ##
  g = g + theme_minimal()
  g = g + theme(panel.border=element_rect(fill=NA,size=1))
  g = g + theme(legend.title=element_blank()) # no legend title
  g = g + theme(axis.title=element_blank()) # no value labels
  g = g + theme(axis.text=element_text(size=16,face='bold',family='Helvetica'))
  g = g + theme(legend.text=element_text(size=14,family='Helvetica'))
  g = g + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
  # g = g + theme(panel.grid.major = element_line(size = 0.25, color = "grey"))
  g
  
  ggsave(outfile)
  
}

plot_circuit1_course(ERKpl300$concentrations.mat, 'plots/circuit1_ERKpl300.pdf', trial=4, plot_legend=T)
plot_circuit1_course(ERKsub300$concentrations.mat, 'plots/circuit1_sub_ERK.pdf', trial=4, plot_legend=T)
