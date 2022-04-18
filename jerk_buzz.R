## R file to plot some figures in the paper 
## "Detection of foraging behavior from accelerometer data using U-Net type convolutional networks"
## Ecological Informatics., 62, 101275, 2021
## [https://www.journals.elsevier.com/ecological-informatics]
## Free version is avaible at [https://arxiv.org/abs/2101.01992]

# Read data ----

library(data.table)
library(arrow)
library(ggplot2)
library(lubridate)


path = 'B:/Codes/Data/Fieldwork 2018'
setwd(path)
source('B:/Codes/MLP_KU/detect_peaks.R')
source('B:/Codes/MLP_KU/njerk.R')
data_path = ''
fs = 100

whale_list = list('Asgeir','Helge18','Kyrri','Nemo','Siggi')
res = list('any', length(whale_list))
threshold_list = list('any', length(whale_list))

dt_info = vector('list', length(whale_list))

for (i in seq.int(length(whale_list))) {
	file_name = paste0('accel-',whale_list[[i]],'.csv.parquet')
	
	dt = read_parquet( paste0(data_path,file_name) )
	setDT(dt)
	dt_info[[i]] = dt[, .(sum(buzz),
	                      {
	                        length(which(diff(.SD$buzz) > 0))
	                      }), by=Dive_no]
}


# Rcpp functions ----

library(Rcpp)

# Choose peaks of X such that distance between 2 consecutive peaks > d [centi-seconds]
cppFunction('
Rcpp::IntegerVector peak_dist(Rcpp::IntegerVector& X, double d, double win_len)
{
	int i, j, k;
	int n = X.size();

	int dist = int(ceil(d/win_len));

	Rcpp::IntegerVector y = Rcpp::IntegerVector(n);


	y[0] = X[0]; k = 1;
	i = 0;
	while (i < n) {
	  for (j = i+1; j<n; j++ ) {
  		if ( (X[j] - X[i]) >= dist ) {
  		  y[k] = X[j];
  		  k++;
  		  break;
  		}
	  }
	  i = j;
	}
	return (y[Rcpp::Range(0, k-1)]);
}
')


# Plot jerk's true positive, precision, and recall (Figure 12) ----

# Function 'jerk_detect' read the data of the narwhal "whale", 
# then try to find all the RMS jerks that is greater than the levels
# Based on that, it detect whether these jerks are true/false positive ,
# then return a dataframe with the true positive, precision, and recall
jerk_detect = function(data_path, whale, break_time = 0.5, t_before = 1, 
											 t_after = 2, steps = 20, levels = NULL, shift_time = 0) 
{

  file_name = paste0('accel-',whale,'.csv.parquet')
  
  dt = read_parquet( paste0(data_path,file_name) )
  setDT(dt)
  dt[,c(6:9):=NULL]
  
  start_buzz = which(diff(dt$buzz) > 0)
  end_buzz = which(diff(dt$buzz) < 0)
    
  if (max(start_buzz%%fs)) {
    start_buzz = start_buzz - start_buzz%%fs
  }
  table(end_buzz%%fs)
  # 0  10  20  30  40  50  60  70  80  90 
  # 96  86  59  68  67  86  72  90  99 113 
  
  x = sqrt(diff(dt$AccX)^2+diff(dt$AccY)^2+diff(dt$AccZ)^2)
  dt$jerk = fs*c(x[1], x)
  
  
  win_len = 20 # centiseconds
  # cut to fit with win_len
  if (nrow(dt) %% win_len) {
    dt = dt[1:(floor(nrow(df)/win_len)*win_len)] }
  
  df = dt[seq(1,nrow(dt),by=win_len), c(1,5:9) ]
  # https://stackoverflow.com/a/30359696
  df$rms_jerk = dt[, sqrt(mean(jerk^2)), by= (1:nrow(dt) - 1) %/% win_len]$V1
  
  
  start_buzz = which(diff(df$buzz) > 0)
  end_buzz = which(diff(df$buzz) < 0)  
  t_s = t_before*fs/win_len
  t_e = t_after*fs/win_len
  
  rms_min = vector('integer', length(start_buzz))
  rms_max = vector('integer', length(start_buzz))
  s = vector('integer', length(start_buzz))
  e = vector('integer', length(end_buzz))
  
  for (i in seq.int(length(start_buzz))) {
    s[i] = start_buzz[i] - t_s
    e[i] = end_buzz[i] + t_e
    rms_min[i] = min(df$rms_jerk[s[i]:e[i]])
    rms_max[i] = max(df$rms_jerk[s[i]:e[i]])
  }
  
  if (shift_time > 0) {
  	s = s + shift_time*fs/win_len
  	e = e + shift_time*fs/win_len
  }
  
  if (is.null(levels)) {
    levels = seq(min(rms_min), max(rms_max), by=steps*fs)
  }

  TP = vector('integer', length(levels))
  FP = vector('integer', length(levels))
  FN = vector('integer', length(levels))
  
  S = ceiling(break_time/(win_len/fs))
  cat( 'No. level:', whale, length(levels), range(levels), '\n' )
  
  for ( i in seq.int( length(levels) ) ) {
    lev = levels[i]
    peaks = which( df$rms_jerk >= lev )
    
    # # Only choose when whales aren't on surface
    # peaks = peaks[which(df[peaks,]$DiveState>0)]
    
    # Only choose when whales are at bottom
    peaks = peaks[which(df[peaks,]$DiveState == 2)]
    peaks = peak_dist(peaks, break_time, win_len/fs)
        
    if ( length(peaks) > 1 ) {
      if ( min(diff(peaks)) < ceiling(break_time/(win_len/fs)) )
        browser()
    }
    
    Buzz = mapply(function(x,y) x:y, s, e)
    buzz_indices = unique(unlist(Buzz))
    # True Positive
    TP[i] = length(intersect(peaks, buzz_indices))
		# False Positive
    FP[i] = length(peaks) - TP[i]
		# False Negative
    if (i == 1) {
    	FN[i] = 0
    } else {
    	FN[i] = TP[1] - TP[i]
    }   

  }
  
  list(
    data.table(whale=rep(whale,length(levels)), level=levels,tp=TP,fp=FP,fn=FN),
    levels
  )
}


for (i in seq.int(length(whale_list))) {
  L = jerk_detect(data_path, whale_list[[i]], break_time = 0.2, 
                  t_before = 0, t_after = 0, 
                  levels = seq(0, 104000, 2000),
                  shift_time = 1.0)
  res[[i]] = L[[1]]
  threshold_list[[i]] = L[[2]]
}


# Plot the curves of all whales
all_whale = data.table(
	whale = rep('All', length(res[[1]]$level) ),
	level = res[[1]]$level,
	tp = rowSums(sapply(res, function(x) x$tp)),
	fp = rowSums(sapply(res, function(x) x$fp)),
	fn = rowSums(sapply(res, function(x) x$fn))
)

all_whale$precision = all_whale$tp/(all_whale$tp+all_whale$fp)
all_whale$recall = all_whale$tp/(all_whale$tp+all_whale$fn)


# Plot the curves for each whale
Res_all = rbindlist(res)
Res_all$precision = Res_all$tp/(Res_all$tp+Res_all$fp)
Res_all$recall = Res_all$tp/(Res_all$tp+Res_all$fn)

Res_all[whale=='Helge18', whale:='Helge']

Res_all = rbind(Res_all, all_whale)

# Change names to (GPS) IDs
name_Id = fread(paste0(data_path,'name_ID.csv'))
Res_all[name_Id, whale := ID , on='whale==Individual']


library(latex2exp)
library(scales)

# https://stackoverflow.com/a/42906139
whale_color = setNames( c(hue_pal()(length(unique(Res_all$whale)) - 1), 'black'),
												unique(Res_all$whale) )
whale_size = setNames( c(rep(1, length(unique(Res_all$whale)) - 1), 2),
											 unique(Res_all$whale) )
theme_set(theme_gray(base_size = 20))

g0 = ggplot(Res_all, aes(x=level,y=tp, size = whale, color=whale)) + 
	scale_color_manual(values = whale_color) +
	scale_size_manual(values = whale_size) + 
	geom_line() + 
  labs(x='Threshold (mG/s)', y='True positive') +
  theme( axis.title.x=element_blank(),
         legend.position="top",
         legend.key.width = unit(30,"mm") ) +
  guides( size = FALSE,
          color = guide_legend(title = '',override.aes = list(size = 2)) ) # + ylim(0,0.2)

g1 = ggplot(Res_all, aes(x=level,y=precision, size = whale, color=whale)) + 
	scale_color_manual(values = whale_color) +
	scale_size_manual(values = whale_size) + 
	geom_line() + 
  labs(x='Threshold (mG/s)', y='Precision') +
  theme( axis.title.x=element_blank(),
         legend.position="none",
         legend.key.width = unit(30,"mm") ) +
  guides( size = FALSE,
          color = guide_legend(title = '',override.aes = list(size = 2)) ) # + ylim(0,0.2)

g2 = ggplot(Res_all, aes(x=level,y=recall,size = whale, color=whale)) + 
	scale_color_manual(values = whale_color) +
	scale_size_manual(values = whale_size) + 
	geom_line() + 
  labs(x='Threshold (mG/s)', y='Recall') +
  theme( legend.position="none",
         legend.key.width = unit(30,"mm") ) +
  guides( size = guide_legend(nrow = 1),
          color = guide_legend(title = '',override.aes = list(size = 2)) )

# Combine all plot together
g0/g1/g2



# Plot dives with jerks (Figure 13) ----
# Plot dive with id "dive_id" from the dataframe of dive data "dt"
# with orange color indicates when jerks happen

win_len = 20
fs = 100
scale_tick = 2

# Choose a dive having buzzes
dt_non_buzz = dt[Dive_no==167]
x = sqrt(diff(dt_non_buzz$AccX)^2+diff(dt_non_buzz$AccY)^2+diff(dt_non_buzz$AccZ)^2)
dt_non_buzz$jerk = fs*c(x[1], x)
# cut to fit with win_len
n = nrow(dt_non_buzz)
dt_non_buzz = dt_non_buzz[1:(floor(n/win_len)*win_len)]
n = nrow(dt_non_buzz)
dt_non_buzz[, RMS_jerk:=rep(sqrt(mean(jerk^2)), each = win_len), by= (1:n - 1) %/% win_len]
dat_non_buzz = dt_non_buzz[seq(1,n,by = win_len)]
dat_non_buzz$Id = 1:nrow(dat_non_buzz)
dat_non_buzz = dat_non_buzz[,c('Id','Depth','RMS_jerk','buzz')]
dat_non_buzz$buzz_type = 'Non-buzzing Dive'

# Choose a dive having no buzzes
dt_buzz = dt[Dive_no==160]
x = sqrt(diff(dt_buzz$AccX)^2+diff(dt_buzz$AccY)^2+diff(dt_buzz$AccZ)^2)
dt_buzz$jerk = fs*c(x[1], x)
# cut to fit with win_len
n = nrow(dt_buzz)
dt_buzz = dt_buzz[1:(floor(n/win_len)*win_len)]
n = nrow(dt_buzz)
dt_buzz[, RMS_jerk:=rep(sqrt(mean(jerk^2)), each = win_len), by= (1:n - 1) %/% win_len]
dat_buzz = dt_buzz[seq(1,n,by = win_len)]
dat_buzz$Id = 1:nrow(dat_buzz)
dat_buzz = dat_buzz[,c('Id','Depth','RMS_jerk','buzz')]
dat_buzz$buzz_type = 'Buzzing Dive'
dat = rbind(dat_non_buzz, dat_buzz)
s = which(diff(dat$buzz) > 0) + 1
e = which(diff(dat$buzz) < 0) + 1

if (length(s) > 0) {
  dat$group = 0
  for (i in seq.int(length(s))) {
    dat[s[i]:e[i],group:=i]
  }
  dt2 = dat[group>0,]
}  


theme_set(theme_gray(base_size = 25))

g1 = ggplot(dat, aes(x = Id) ) + 
  geom_line( aes(y = Depth, color=factor(buzz), group=1), 
             size=1, show.legend=F ) +
  scale_x_continuous(labels = function(x) x*win_len/(60*fs), 
                     breaks = function(x) seq(0, floor(x[2]), 
                                              by = scale_tick*60*fs/win_len) ) +
  scale_y_continuous(trans = 'reverse') +
  scale_color_manual(values = c("0" = "black", "1" = "#D55E00") ) +
  scale_size_manual(values = c(0.5, 2)) +
  labs( y = "Depth (m)", x = "Time (minutes)" ) + 
  geom_ribbon(data = dt2 , fill = "#D55E00",
              aes(ymax = max(Depth)+25, ymin = Depth, group = group)) +
  facet_wrap( ~ buzz_type, ncol=1) + theme(axis.title.x=element_blank() )


g2 = ggplot(dat, aes(x = Id) ) + 
  geom_line( aes(y = RMS_jerk, color=factor(buzz), group=1), 
             size=1, show.legend=F ) +
  scale_x_continuous(labels = function(x) x*win_len/(60*fs), 
                     breaks = function(x) seq(0, floor(x[2]), 
                                              by = scale_tick*60*fs/win_len) ) +
  scale_color_manual(values = c("0" = "black", "1" = "#D55E00") ) +
  scale_size_manual(values = c(0.5, 2)) +
  labs( y = "RMS jerk (mG/s)", x = "Time (minutes)" ) +
  geom_ribbon(data = dt2 , fill = "#D55E00",
              aes(ymin = -25, ymax = RMS_jerk, group = group)) +
  facet_wrap( ~ buzz_type, ncol=1) + theme(axis.title.x=element_blank() )


## Combine plot together
gp = g1|g2 
gp + plot_annotation(  caption = 'Time (minutes)', 
                       theme = theme(plot.caption = element_text(size = 25, hjust = 0.55) ) )
