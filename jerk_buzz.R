# ----

library(data.table)
library(arrow)
library(ggplot2)

machine = 2

if (machine == 1) {
  path = 'D:/codes/MLP_ku'
  setwd(path)
  source('detect_peaks.R')
  source('njerk.R') 
  data_path = 'data/'
} else if (machine == 2) {
  path = 'B:/Codes/Data/Fieldwork 2018'
  setwd(path)
  source('B:/Codes/MLP_KU/detect_peaks.R')
  source('B:/Codes/MLP_KU/njerk.R')
  data_path = ''
}

fs = 100

whale_list = list('Asgeir','Helge18','Kyrri','Nemo','Siggi')
res = list('any', length(whale_list))
threshold_list = list('any', length(whale_list))

dt_info = vector('list', length(whale_list))

for (i in seq.int(length(whale_list))) {
	file_name = paste0('accel-',whale_list[[i]],'.csv.parquet')
	
	dt = read_parquet( paste0(data_path,file_name) )
	setDT(dt)
	
	# dt_info[[i]] = dt[Dive_no > 0, sum(buzz), by=Dive_no]
	# 
	# [https://stackoverflow.com/a/2589834]9
	dt_info[[i]] = dt[, .(sum(buzz),
	                      {
	                        length(which(diff(.SD$buzz) > 0))
	                      }), by=Dive_no]
}


range(dt$AccX)
range(dt$AccY)
range(dt$AccZ)


library(lubridate)

seconds_to_period(sum(sapply( dt_info, function(x) sum(x$V1) )))
sum(sapply( dt_info, function(x) sum(x$V2) ))


sum(sapply(dt_info, nrow))
sum(sapply(dt_info, function(x) nrow(x[V1 > 0,]) ))
sum(sapply(dt_info, function(x) nrow(x[V1 > 0,]) ))/sum(sapply(dt_info, nrow))



# whale = whale_list[[1]]


debug_jerk_detect = function(data_path, whale, t_before = 1, t_after = 2) 
{
  
  file_name = paste0('accel-',whale,'.csv.parquet')
  
  dt = read_parquet( paste0(data_path,file_name) )
  setDT(dt)
  dt[,c(6:9):=NULL]
  
  start_buzz = which(diff(dt$buzz) > 0)
  end_buzz = which(diff(dt$buzz) < 0)
  
  end_buzz - start_buzz
  
  if (max(start_buzz%%fs)) {
    start_buzz = start_buzz - start_buzz%%fs
  }
  table(end_buzz%%fs)
  # 0  10  20  30  40  50  60  70  80  90 
  # 96  86  59  68  67  86  72  90  99 113 
  
  
  dt$jerk = njerk(matrix(c(dt$AccX,dt$AccY,dt$AccZ), 
                         byrow = T, ncol = 3), 
                  sampling_rate = fs)/fs
  
  
  win_len = 20
  
  # cut to fit with win_len
  if (nrow(dt) %% win_len) {
    dt = dt[1:(floor(nrow(df)/win_len)*win_len)] }
  
  df = dt[seq(1,nrow(dt),by=win_len), c(1,5:9) ]
  # https://stackoverflow.com/a/30359696
  df$rms_jerk = dt[, sqrt(mean(jerk^2)), by= (1:nrow(dt) - 1) %/% win_len]$V1
  
  gc()
  
  
  start_buzz = which(diff(df$buzz) > 0)
  end_buzz = which(diff(df$buzz) < 0)
  
  (end_buzz - start_buzz)*win_len
  
  break_time = 0.5
  
  t_s = t_before*fs/win_len
  t_e = t_after*fs/win_len
  
  
  rms_min = vector('integer', length(start_buzz))
  rms_max = vector('integer', length(start_buzz))
  s = vector('integer', length(start_buzz))
  e = vector('integer', length(end_buzz))
  
  
  length(start_buzz)
  length(end_buzz)
  
  for (i in seq.int(length(start_buzz))) {
    s[i] = start_buzz[i] - t_s
    e[i] = end_buzz[i] + t_e
    rms_min[i] = min(df$rms_jerk[s[i]:e[i]])
    rms_max[i] = max(df$rms_jerk[s[i]:e[i]])
    
  }
  
  levels = seq(min(rms_min), max(rms_max), by=20)
  
  TP = vector('integer', length(levels))
  FP = vector('integer', length(levels))
  
  i = 1
  Y = vector('list', length(s))
  Z = vector('list', length(s))
  for (i in seq.int(length(levels))) {
    lev = levels[i]
    peaks = which( df$rms_jerk > lev )
    # Only choose when whales aren't on surface
    peaks = peaks[which(df[peaks,]$Dive_no>0)]
    
    # for (k in seq.int(length(s))) {
    #   # X = peaks[which((peaks - s[k] ) >= 0 & (peaks - e[k] ) <= 0)]
    #   Y[[k]] = intersect(peaks, as.vector(mapply(function(x,y) x:y, s[k], e[k])))
    #   # Y[[k]] = as.vector(mapply(function(x,y) x:y, s[k], e[k]))
    #   # if (!identical(X,Y[[k]])) {
    #   #   browser()
    #   # }
    # }

    q = mapply(function(x,y) x:y, s, e)

    for (k in seq.int(length(s))) {
      Z[[k]] = intersect(peaks, q[[k]])
      # if (!identical(Y[[k]], Z[[k]])) {
      #   browser()
      # }
    }
    w0 = unlist(Z)
    
    buzz_indices = unlist(q)
    w = intersect(peaks, buzz_indices)

    p = which(w0[1:length(w)]!=w)[1]
    
    pos = which( p < cumsum(sapply(q, length)) )[1]
    
    # 17 & 18
    
    
      browser()
  }
  
  data.table(whale=rep(whale,length(levels)), 
             level=levels,tp=TP,fp=FP)
  
}

# Rcpp functions ----

library(Rcpp)

# [http://adv-r.had.co.nz/Rcpp.html]
# [https://stat.ethz.ch/pipermail/r-help/2014-October/422684.html]
# [https://teuder.github.io/rcpp4everyone_en/080_vector.html]

# Choose peaks s.t. distance between 2 consecutive peaks > d [centi-seconds]
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

  // [https://stackoverflow.com/a/47246770]
	return (y[Rcpp::Range(0, k-1)]);
	
}
')



cppFunction('
int intersect_cpp(Rcpp::IntegerVector& x, Rcpp::IntegerVector s, Rcpp::IntegerVector e)
{
	int i, cnt = 0, k;
	int n = x.size();
	int m = s.size();

	// Rcpp::IntegerVector y = Rcpp::IntegerVector(n);

	k = 0;
  for (i = 0; i < n; i++ ) {
		if ( (x[i] >= s[k]) && (x[i] <= e[k]) ) {
		  cnt++;
		} else {
		  k++ ;
		  if (k>=m) 
		    break;
		}
  }

	return cnt ;
	
}
')

# Plot jerk's precision & recall ----

jerk_detect = function(data_path, whale, break_time = 0.5, t_before = 1, 
											 t_after = 2, steps = 20, levels = NULL, shift_time = 0) 
{

  file_name = paste0('accel-',whale,'.csv.parquet')
  
  dt = read_parquet( paste0(data_path,file_name) )
  setDT(dt)
  dt[,c(6:9):=NULL]
  
  start_buzz = which(diff(dt$buzz) > 0)
  end_buzz = which(diff(dt$buzz) < 0)
  
  end_buzz - start_buzz
  
  if (max(start_buzz%%fs)) {
    start_buzz = start_buzz - start_buzz%%fs
  }
  table(end_buzz%%fs)
  # 0  10  20  30  40  50  60  70  80  90 
  # 96  86  59  68  67  86  72  90  99 113 
  
  # dt$jerk = njerk(matrix(c(dt$AccX,dt$AccY,dt$AccZ), 
  #                        byrow = T, ncol = 3), sampling_rate = fs)
  x = sqrt(diff(dt$AccX)^2+diff(dt$AccY)^2+diff(dt$AccZ)^2)
  dt$jerk = fs*c(x[1], x)
  
  
  win_len = 20 # centiseconds

  # cut to fit with win_len
  if (nrow(dt) %% win_len) {
    dt = dt[1:(floor(nrow(df)/win_len)*win_len)] }
  
  df = dt[seq(1,nrow(dt),by=win_len), c(1,5:9) ]
  # https://stackoverflow.com/a/30359696
  df$rms_jerk = dt[, sqrt(mean(jerk^2)), by= (1:nrow(dt) - 1) %/% win_len]$V1
  
  # gc()
  
  start_buzz = which(diff(df$buzz) > 0)
  end_buzz = which(diff(df$buzz) < 0)
  (end_buzz - start_buzz)*win_len
  
  t_s = t_before*fs/win_len
  t_e = t_after*fs/win_len
  
  rms_min = vector('integer', length(start_buzz))
  rms_max = vector('integer', length(start_buzz))

  s = vector('integer', length(start_buzz))
  e = vector('integer', length(end_buzz))

  length(start_buzz)
  length(end_buzz)
  
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
  
  # if (df$whale[1] == 'Asgeir') {
  #   browser()
  # }
    
  for ( i in seq.int( length(levels) ) ) {
    lev = levels[i]
    peaks = which( df$rms_jerk >= lev )
    
    # # Only choose when whales aren't on surface
    # peaks = peaks[which(df[peaks,]$DiveState>0)]
    
    # Only choose when whales are at bottom
    peaks = peaks[which(df[peaks,]$DiveState == 2)]
    peaks = peak_dist(peaks, break_time, win_len/fs)
    
    # if ( is.infinite(min(diff(peaks))) ) {
    #   browser() 
    # }
    
    if ( length(peaks) > 1 ) {
      if ( min(diff(peaks)) < ceiling(break_time/(win_len/fs)) )
        browser()
    }
    
    # TP[i] = intersect_cpp(peaks, s, e)
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
    
    # browser()
    # split(buzz_indices, cumsum(c(1, diff(buzz_indices) != 1)))
    # FN[i] = sum( ceiling( sapply(Buzz, length)/S ) ) - TP[i]

  }
  
  list(
    data.table(whale=rep(whale,length(levels)), level=levels,tp=TP,fp=FP,fn=FN),
    levels
  )
}


for (i in seq.int(length(whale_list))) {
  # res[[i]] = debug_jerk_detect(data_path, whale_list[[i]])
  L = jerk_detect(data_path, whale_list[[i]], break_time = 0.2, 
                  t_before = 0, t_after = 0, 
                  levels = seq(0, 104000, 2000),
                  shift_time = 1.0)
  res[[i]] = L[[1]]
  threshold_list[[i]] = L[[2]]
}


# # Threshold plot
# theme_set(theme_gray(base_size = 30))
# d = data.frame(x = unlist(threshold_list), 
#                 grp = rep(unlist(whale_list),
#                           times = sapply(threshold_list,length)))
# d[grp=='Helge18', grp:='Helge']
# ggplot(d,aes(x = grp, y = x)) + geom_boxplot() + labs( x = 'Whale', y = 'Threshold (mG/s)')


# Plot all whales
all_whale = data.table(
	whale = rep('All', length(res[[1]]$level) ),
	level = res[[1]]$level,
	tp = rowSums(sapply(res, function(x) x$tp)),
	fp = rowSums(sapply(res, function(x) x$fp)),
	fn = rowSums(sapply(res, function(x) x$fn))
)

all_whale$precision = all_whale$tp/(all_whale$tp+all_whale$fp)
all_whale$recall = all_whale$tp/(all_whale$tp+all_whale$fn)



# Plot each whale
Res_all = rbindlist(res)
# Res_all$whale = unlist(Res_all$whale)
# Res_all$precision = log10(Res_all$tp/(Res_all$tp+Res_all$fp))
Res_all$precision = Res_all$tp/(Res_all$tp+Res_all$fp)
Res_all$recall = Res_all$tp/(Res_all$tp+Res_all$fn)

Res_all[whale=='Helge18', whale:='Helge']

Res_all = rbind(Res_all, all_whale)

# Change names to (GPS) IDs
name_Id = fread(paste0(data_path,'name_ID.csv'))
Res_all[name_Id, whale := ID , on='whale==Individual']


library(ggplot2)
library(latex2exp)
library(patchwork)

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
  # labs(x=TeX('Level $\\left(m\\cdot s^{-2} \\right)$'),y='Precision') +
  labs(x='Threshold (mG/s)', y='True positive') +
  # theme( legend.position="top" ) +
  # guides(size = guide_legend(nrow = 1), color=guide_legend(""))
  theme( axis.title.x=element_blank(),
         legend.position="top",
         legend.key.width = unit(30,"mm") ) +
  guides( size = FALSE,
          color = guide_legend(title = '',override.aes = list(size = 2)) ) # + ylim(0,0.2)

g1 = ggplot(Res_all, aes(x=level,y=precision, size = whale, color=whale)) + 
	scale_color_manual(values = whale_color) +
	scale_size_manual(values = whale_size) + 
	geom_line() + 
  # labs(x=TeX('Level $\\left(m\\cdot s^{-2} \\right)$'),y='Precision') +
  labs(x='Threshold (mG/s)', y='Precision') +
  # theme( legend.position="top" ) +
  # guides(size = guide_legend(nrow = 1), color=guide_legend(""))
  theme( axis.title.x=element_blank(),
         # legend.position="top",
         legend.position="none",
         legend.key.width = unit(30,"mm") ) +
  guides( size = FALSE,
          color = guide_legend(title = '',override.aes = list(size = 2)) ) # + ylim(0,0.2)


g2 = ggplot(Res_all, aes(x=level,y=recall,size = whale, color=whale)) + 
	scale_color_manual(values = whale_color) +
	scale_size_manual(values = whale_size) + 
	geom_line() + 
  # labs(x=TeX('Level $\\left(mG\\cdot s^{-1} \\right)$'),y='Recall') +
  labs(x='Threshold (mG/s)', y='Recall') +
  theme( legend.position="none",
         legend.key.width = unit(30,"mm") ) +
  guides( size = guide_legend(nrow = 1),
          color = guide_legend(title = '',override.aes = list(size = 2)) )

  
# g1/g2
g0/g1/g2



# range(df$rms_jerk[s:e])
# # [1]  38.46218 278.24732
# 
# data.frame(detect_peaks(df$rms_jerk[s:e], fs/win_len, 
#              bktime = break_time, 
#              thresh = 44.0, plot_peaks=F)[1:4])



#

# small dataset ----

whale = whale_list[[1]]
file_name = paste0('accel-',whale,'.csv.parquet')
dt = read_parquet( paste0(data_path,file_name) )
setDT(dt)

ID = which(diff(dt$buzz) > 0)
id = ID[2:6]
m = length(id)
df = dt[(id[1]-100):(id[m]+299) ,]


buzz_dive = unique(dt[ID+1,Dive_no])
buzz_dive
#  81  86  87  89  91  92  93  98 107 117 118 119 120 121 122 126 127 128 130 131 133 
# 134 135 136 139 145 147 148 149 157 158 159
# 160 165 166 169 170 171 172 173 174 175 176 178 179 180 181 182 183 184 185 186 193


no_buzz_dive = setdiff(seq.int(buzz_dive[1], tail(buzz_dive,1)), buzz_dive)

# for (i in no_buzz_dive) {
#   print(max(dt[Dive_no == i,Depth]))
# }
  

no_buzz_dive_depth = dt[Dive_no %in% no_buzz_dive, max(Depth), by=Dive_no ]
names(no_buzz_dive_depth) = c('Dive_no','Max_depth')

# View(no_buzz_dive_depth)

buzz_dive_depth = dt[Dive_no %in% buzz_dive, max(Depth), by=Dive_no ]
names(buzz_dive_depth) = c('Dive_no','Max_depth')

# View(buzz_dive_depth)

# ----

i = min( unique(abs(dt$Dive_no)) ) + 1

while(T) {
  df = dt[Dive_no==i]
  if (max(df$buzz) > 0) {
    break
  }

  i = i+1
}


# Plot dives with jerks ----

library(ggplot2)
library(patchwork)

theme_set(theme_gray(base_size = 22))


plot_dive = function(dt, dive_id, win_len = 20, scale_tick=2) 
{
  
  df = dt[Dive_no==dive_id]

  fs = 100
  
  # fcut = 0.25
  #
  # source('B:/Codes/MLP_KU/fir_nodelay.R')
  # highpass_freq = fcut/(fs/2)
  #
  # https://dsp.stackexchange.com/a/3183
  # https://stackoverflow.com/a/63733511
  # 
  # Ax_filt = fir_nodelay(x=dt$AccX, n=51, qual='high', fc = highpass_freq )
  # Ay_filt = fir_nodelay(x=dt$AccY, n=51, qual='high', fc = highpass_freq )
  # Az_filt = fir_nodelay(x=dt$AccZ, n=51, qual='high', fc = highpass_freq )
  # Ax_filt = fir_nodelay(x=df$AccX, n=101, qual='high', fc = highpass_freq )
  # Ay_filt = fir_nodelay(x=df$AccY, n=101, qual='high', fc = highpass_freq )
  # Az_filt = fir_nodelay(x=df$AccZ, n=101, qual='high', fc = highpass_freq )
  # Amatrux = matrix(c(Ax_filt,Ay_filt,Az_filt), byrow = TRUE, ncol = 3)
  
  Amatrux = matrix(c(df$AccX,df$AccY,df$AccZ), byrow = TRUE, ncol = 3)
  head(Amatrux)
  
  J = njerk(Amatrux, sampling_rate = fs)
  
  df$jerk = J
  
  
  # https://stackoverflow.com/a/30359696
  
  # cut to fit with win_len
  df = df[1:(floor(nrow(df)/win_len)*win_len)]
  
  df[, RMS_jerk:=rep(sqrt(mean(jerk^2)), each = win_len), 
        by= (1:nrow(df) - 1) %/% win_len]

  dat = df[seq(1,nrow(df),by = win_len)]
  
  # Random grid search
  # has a probability of 95% of finding a combination of parameters within the 5% optima with only 60 iterations
  # https://stats.stackexchange.com/a/209409
  # thresh = 80000
  # peaks = detect_peaks(J, fs/win_len, bktime = 5, thresh = thresh, plot_peaks=F)
  
  # pks1 = detect_peaks(J, fs, bktime = 5, thresh = 110000, plot_peaks=T)

  if (max(dat$buzz)) {
    plot_title = 'Buzzing dive'
  } else {
    plot_title = 'Non-buzzing dive'
  }
  
  dat$Id = 1:nrow(dat)
  
  g1 = ggplot(dat, aes(x = Id) ) + 
    geom_line( aes(y = RMS_jerk, color=factor(buzz),group=1), 
               size=1, show.legend = F ) + 
    scale_color_manual(values = c("0" = "black", "1" = "#D55E00") ) +
    scale_x_continuous(labels = function(x) x*win_len/(60*fs), # https://stackoverflow.com/a/54678868
                       breaks = function(x) seq(0, floor(x[2]), 
                                                by = scale_tick*60*fs/win_len)) +
    scale_size_manual(values = c(0.5, 1)) +
    labs( y = "RMS jerk", x = "" )
  
  
  g2 = ggplot(dat, aes(x = Id) ) + 
    geom_line( aes(y = Depth, color=factor(buzz), group=1), 
               size=1, show.legend=F ) +
    scale_x_continuous(labels = function(x) x*win_len/(60*fs), 
                       breaks = function(x) seq(0, floor(x[2]), 
                                                by = scale_tick*60*fs/win_len) ) +
    scale_y_continuous(trans = 'reverse') +
    scale_color_manual(values = c("0" = "black", "1" = "#D55E00") ) +
    scale_size_manual(values = c(0.5, 2)) +
    labs( y = "Depth (m)", x = "Time (minutes)" )
  
  
  s = which(diff(dat$buzz) > 0) + 1
  e = which(diff(dat$buzz) < 0) + 1
  
  if (length(s) > 0) {
    dat$group = 0
    for (i in seq.int(length(s))) {
      dat[s[i]:e[i],group:=i]
    }
    dt2 = dat[group>0,]
    
    g1 = g1 + geom_ribbon(data = dt2 , 
                aes(ymin = -25, ymax = RMS_jerk, group = group),
                fill = "#D55E00")
    
    g2 = g2 + geom_ribbon(data = dt2 , 
                aes(ymax = max(Depth)+25, ymin = Depth, group = group),
                fill = "#D55E00")
  }
  

  # https://stackoverflow.com/a/34805932
  g1 = g1 + ggtitle(plot_title) + 
    theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5 ))
  

  g2 = g2 + ggtitle(plot_title) + 
    theme(plot.title = element_text(margin = margin(t = 10, b = -20), hjust = 0.5 ))
  
  
  g1/g2
}


# plot_dive(dt, 86)
# plot_dive(dt, 87)
# plot_dive(dt, 89)

  p1 = plot_dive(dt, 167)
  plot_dive(dt, 162)
  plot_dive(dt, 109)
  plot_dive(dt, 163)


plot_dive(dt, 165)
plot_dive(dt, 159)
p2 = plot_dive(dt, 160)
plot_dive(dt, 170)


p1/p2[[1]]/p2[[2]]



# https://rdrr.io/cran/cvTools/man/cvTuning.html

# Facet wrap ----

win_len = 20
fs = 100
scale_tick = 2

dt_non_buzz = dt[Dive_no==167]
# Amatrux = matrix(c(dt_non_buzz$AccX,dt_non_buzz$AccY,dt_non_buzz$AccZ), 
#                  byrow = TRUE, ncol = 3)
# dt_non_buzz$jerk = njerk(Amatrux, sampling_rate = fs)
x = sqrt(diff(dt_non_buzz$AccX)^2+diff(dt_non_buzz$AccY)^2+diff(dt_non_buzz$AccZ)^2)
dt_non_buzz$jerk = fs*c(x[1], x)
# cut to fit with win_len
n = nrow(dt_non_buzz)
dt_non_buzz = dt_non_buzz[1:(floor(n/win_len)*win_len)]
n = nrow(dt_non_buzz)
dt_non_buzz[, RMS_jerk:=rep(sqrt(mean(jerk^2)), each = win_len), by= (1:n - 1) %/% win_len]
dat_non_buzz = dt_non_buzz[seq(1,n,by = win_len)]
dat_non_buzz$Id = 1:nrow(dat_non_buzz)
  # dat_non_buzz$Depth = -dat_non_buzz$Depth
dat_non_buzz = dat_non_buzz[,c('Id','Depth','RMS_jerk','buzz')]
  # dat_non_buzz = melt(dat_non_buzz, id.vars = c("Id",'buzz'))
dat_non_buzz$buzz_type = 'Non-buzzing Dive'


dt_buzz = dt[Dive_no==160]
# Amatrux = matrix(c(dt_buzz$AccX,dt_buzz$AccY,dt_buzz$AccZ), 
#                  byrow = TRUE, ncol = 3)
# dt_buzz$jerk = njerk(Amatrux, sampling_rate = fs)
x = sqrt(diff(dt_buzz$AccX)^2+diff(dt_buzz$AccY)^2+diff(dt_buzz$AccZ)^2)
dt_buzz$jerk = fs*c(x[1], x)
# cut to fit with win_len
n = nrow(dt_buzz)
dt_buzz = dt_buzz[1:(floor(n/win_len)*win_len)]
n = nrow(dt_buzz)
dt_buzz[, RMS_jerk:=rep(sqrt(mean(jerk^2)), each = win_len), by= (1:n - 1) %/% win_len]
dat_buzz = dt_buzz[seq(1,n,by = win_len)]
dat_buzz$Id = 1:nrow(dat_buzz)
  # dat_buzz$Depth = -dat_buzz$Depth
dat_buzz = dat_buzz[,c('Id','Depth','RMS_jerk','buzz')]
  # dat_buzz = melt(dat_buzz, id.vars = c("Id",'buzz'))
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


# g1 = g1 + geom_ribbon(data = dt2 , 
#                       aes(ymin = -25, ymax = RMS_jerk, group = group),
#                       fill = "#D55E00")
# 
# g2 = g2 + geom_ribbon(data = dt2 , 
#                       aes(ymax = max(Depth)+25, ymin = Depth, group = group),
#                       fill = "#D55E00")

theme_set(theme_gray(base_size = 25))

# ggplot(dat, aes(x = Id) ) + 
#   geom_line( aes(y = value, color=factor(buzz), group=1), 
#              size=1, show.legend=F ) +
#   scale_x_continuous(labels = function(x) x*win_len/(60*fs), 
#                      breaks = function(x) seq(0, floor(x[2]), 
#                                               by = scale_tick*60*fs/win_len) ) +
#   # scale_y_continuous(labels = function(x) -x) +
#   scale_color_manual(values = c("0" = "black", "1" = "#D55E00") ) +
#   scale_size_manual(values = c(0.5, 2)) +
#   labs( y = "Depth (m)", x = "Time (minutes)" ) +
#   geom_ribbon(data = dt2 , fill = "#D55E00", 
#               # aes(ymax = max(Depth)+25, ymin = Depth, group = group)) +
#               aes(ymin = -25, ymax = value, group = group)) +
#   facet_wrap( variable ~ buzz_type, scales = "free_y", 
#              nrow = 2)


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


library(patchwork)

gp = g1|g2 
gp + plot_annotation(  caption = 'Time (minutes)', 
                       theme = theme(plot.caption = element_text(size = 25, hjust = 0.55) ) )

