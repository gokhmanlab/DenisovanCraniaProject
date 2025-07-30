compare_specimen_percentiles <- function(data, specimen, predictions,
                                         normalize = FALSE,
                                         norm_method = 'pensize',
                                         maxNApercent = 100,
                                         calc.pval = 'U-test',
                                         vs = 'both',
                                         percentile.dist,
                                         remove.outliers = TRUE,
                                         cor.adj = TRUE,
                                         adj.m = c()) {
  
  
  get.per <- function(m, specimen, group, remove.outliers = remove.outliers) {
    if (!(specimen %in% m$specimen)) {
      cat('\t\tN/A\n')
      return(NA)
    }
    test.value <- m[m$specimen == specimen, 'value']
    cat(paste('\t\tvalue:', round(test.value, 3)), '\n')
    
    m <- m[m$specimen != specimen, ]
    compare.values <- if (group == 'AMH') {
      m[m$group %in% c('UPS', 'EHS'), 'value']
    } else {
      m[m$group == group, 'value']
    }
    
    if (sum(!is.na(compare.values)) < 5) {
      cat("Skipping iteration due to insufficient non-NA values\n")
      return(NA)
    }
    
    if (remove.outliers) {
      compare.values <- compare.values[!compare.values %in% boxplot.stats(compare.values)$out]
    }
    
    if (percentile.dist == 'quantile') {
      compare.values <- sort(compare.values)
      N <- length(compare.values)
      Quant <- (1:N - 0.5) / N
      group.cdf <- approxfun(compare.values, Quant, yleft = 0, yright = 1, ties = 'mean')
      return(100 * group.cdf(test.value))
    } else if (percentile.dist == 'normal') {
      SD <- sd(compare.values, na.rm = TRUE)
      MEAN <- mean(compare.values, na.rm = TRUE)
      Z <- Normal(0, 1)
      normal.value <- (test.value - MEAN) / SD
      c.prob <- round(cdf(Z, normal.value), 3)
      cat('\t\tPercentile:', c.prob * 100, '\n')
      return(c.prob * 100)
    }
  }
  
  cat(paste('\n', specimen, '\n'))
  groups <- c('AMH', 'NE')
  
  # Confirm that specimen.group is uniquely defined
  specimen.group <- unique(data$group[data$specimen == specimen])
  if (length(specimen.group) != 1) stop("specimen.group not uniquely defined")
  
  relevant.measurements <- predictions %>%
    filter(include.in.analysis == 'yes' & loc %in% c(1, 2)) %>%
    pull(long.name)
  
  percentiles <- data.frame(matrix(ncol = 6, nrow = 0))
  names(percentiles) <- c('related.prediction', 'prediction', 'Vs', 'direction', 'percentile', 'loc')
  
  measurement_data_list <- list()
  
  for (measurement in relevant.measurements) {
    fixed.measurement <- measurement
    cat(paste('\t', fixed.measurement, '\n'))
    
    measurement_clean <- gsub('[() :,/\\-]', '.', measurement)
    
    m <- get.m(data, fixed.measurement, predictions,
               normalize = normalize,
               method = norm_method,
               TS = specimen.group,
               adj.m = adj.m,
               cor.adj = cor.adj)
    
    measurement_data_list[[measurement_clean]] <- m
    
    prediction_row <- predictions[predictions$long.name == fixed.measurement, ]
    
    for (group in groups) {
      percentile <- tryCatch({
        get.per(m, specimen, group, remove.outliers = remove.outliers)
      }, error = function(e) {
        browser()
      })
      
      percentiles[nrow(percentiles) + 1, ] <- c(
        related.prediction = as.character(prediction_row$related.prediction),
        prediction = fixed.measurement,
        Vs = group,
        direction = as.character(prediction_row[[group]]),
        percentile = percentile,
        loc = prediction_row$loc
      )
    }
  }
  
  dir.loc <- hash()
  dir.loc[['less']] <- 1
  dir.loc[['undetected']] <- 2
  dir.loc[['greater']] <- 3
  percentiles$direction <- values(dir.loc, keys = percentiles$direction)
  
  percentiles <- percentiles[!is.na(percentiles$percentile), ]
  percentiles$percentile <- as.numeric(percentiles$percentile)
  percentiles$prediction <- NULL
  
  aggregated_percentiles <- aggregate(percentile ~ ., percentiles, mean)
  aggregated_percentiles <- subset(aggregated_percentiles, direction != 2)
  
  if (vs == 'AMH') {
    aggregated_percentiles <- subset(aggregated_percentiles, Vs == 'AMH')
  } else if (vs == 'NE') {
    aggregated_percentiles <- subset(aggregated_percentiles, Vs == 'NE')
  } else if (vs == 'both.strict') {
    aggregated_percentiles <- aggregated_percentiles[!duplicated(aggregated_percentiles[, c(1, 3)]), ]
  }
  
  if (!is.na(calc.pval)) {
    aggregated_percentiles$percentile <- aggregated_percentiles$percentile / 100
    
    directional_deviation <- numeric(nrow(aggregated_percentiles))
    for (i in seq_len(nrow(aggregated_percentiles))) {
      directional_deviation[i] <- if (aggregated_percentiles[i, 'direction'] == 3) {
        aggregated_percentiles[i, 'percentile'] - 0.5
      } else {
        (1 - aggregated_percentiles[i, 'percentile']) - 0.5
      }
    }
    aggregated_percentiles$d <- directional_deviation
    
    n <- nrow(aggregated_percentiles)
    mean_d <- round(mean(directional_deviation), 3)
    n_matched <- sum(directional_deviation > 0)
    sum_abs_d <- sum(abs(directional_deviation))
    EKB <- data[data$specimen == specimen, 'EKB.vis']
    cranial_capacity <- data[data$specimen == specimen, 'Cranial.capacit']
    
    binom.p.val <- wilcox.p.val <- NA
    
    if (calc.pval == 'distance mean') {
      z <- sum(directional_deviation) / sqrt((1 / 12) / n)
      p.val <- 1 - pnorm(z)
    } else if (calc.pval == 'U-test') {
      wilcox.p.val <- round(wilcox.test(directional_deviation, alternative = 'g', mu = 0)$p.value, 7)
    } else if (calc.pval == 'binomial') {
      binom.p.val <- round(binom.test(sum(directional_deviation > 0), n,
                                      p = 0.5, alternative = 'g')$p.value, 7)
      wilcox.p.val <- round(wilcox.test(directional_deviation, alternative = 'g', mu = 0)$p.value, 7)
    } else if (calc.pval == 'mid.range.binomial') {
      x <- sum(directional_deviation > 0)
      if (x == n) {
        binom.p.val <- 0.5 * round(binom.test(x, n, p = 0.5, alternative = 'g')$p.value, 7)
      } else {
        p1 <- round(binom.test(x, n, p = 0.5, alternative = 'g')$p.value, 7)
        p2 <- round(binom.test(x + 1, n, p = 0.5, alternative = 'g')$p.value, 7)
        binom.p.val <- 0.5 * (p1 - p2) + p2
      }
      
      wilcox.test <- wilcoxsign_test(
        directional_deviation ~ rep(0, length(directional_deviation)),
        distribution = "exact", alternative = 'g'
      )
      wilcox.p.val <- round(wilcox.test@distribution@pvalue(wilcox.test@statistic@teststatistic), 7)
    }
    
    pval <- data.frame(
      specimen = specimen,
      group = as.character(specimen.group),
      binom.p.val = binom.p.val,
      binom.score = -log10(binom.p.val),
      wilcox.p.val = wilcox.p.val,
      wilcox.score = -log10(wilcox.p.val),
      comb_score = sqrt((-log10(binom.p.val))^2 + (-log10(wilcox.p.val))^2),
      mean.d = mean_d,
      sum.abs.d = sum_abs_d,
      total_predictions = n,
      match_d_over_0.5 = n_matched,
      match_rate = round(n_matched / n, 3),
      perm_success_binom = NA,
      perm_success_matched_d = NA,
      perm_success_U = NA,
      perm_success_mean_d = NA,
      perm_success_dist = NA,
      permutations = NA,
      normalize = normalize,
      norm_method = norm_method,
      EKB = EKB,
      cranial_capacity = cranial_capacity,
      greater_than_AMH = sum(aggregated_percentiles$Vs == 'AMH' & aggregated_percentiles$direction == 3),
      less_than_AMH = sum(aggregated_percentiles$Vs == 'AMH' & aggregated_percentiles$direction == 1),
      greater_than_NE = sum(aggregated_percentiles$Vs == 'NE' & aggregated_percentiles$direction == 3),
      less_than_NE = sum(aggregated_percentiles$Vs == 'NE' & aggregated_percentiles$direction == 1)
    )
    
    return(list(pval, measurement_data_list, aggregated_percentiles))
  }
}


run_permutations_for_specimen <- function(data, specimen, predictions, measurements,
                                          maxNApercent = 100, normalize = FALSE,
                                          calc.pval = 'U-test', norm_method = 'pensize',
                                          remove.outliers = TRUE, percentile.dist = 'quantile',
                                          cor.adj = TRUE, adj.m = c(),
                                          permutation.number = 1000) {
  
  # Filter specimen measurements
  spec <- filter(data, specimen == !!specimen)
  specimen.measurements <- measurements[measurements %in% colnames(spec)[!is.na(spec)]]
  
  # Keep only predictions available in this specimen
  perm.predictions <- predictions %>%
    filter(long.name %in% specimen.measurements)
  
  # Prepare permutations results table
  permutations <- as.data.frame(matrix(nrow = permutation.number, ncol = 28))
  colnames(permutations) <- c(
    'specimen', 'group', 'binom.p.val', 'binom.score', 'wilcox.p.val',
    'wilcox.score', 'comb_score', 'mean.d', 'sum.abs.d',
    'total_predictions', 'match_d_over_0.5', 'match_rate',
    'perm_success_binom', 'perm_success_matched_d', 'perm_success_U',
    'perm_success_mean_d', 'perm_success_dist', 'permutations',
     'normalize', 'norm_method', 'EKB',
    'cranial_capacity',
    'greater_than_AMH', 'less_than_AMH',
    'greater_than_NE', 'less_than_NE',
    'sampled_predictions','available.measurements'
  )
  relevant.perm.count <- 0
  per.sum.abs.d <- numeric(0)
  
  while (relevant.perm.count < permutation.number) {
    cat('\n', relevant.perm.count)
    
    #perm.predictions$long.name <- sample(specimen.measurements, size = nrow(perm.predictions))
    
    #### Generate predictions that account for angle/distance in the original predictions
    # 1. Split based on type
    perm.dist <- perm.predictions %>% filter(type == "distance")
    perm.angle <- perm.predictions %>% filter(type == "angle")
    
    # 2. Subset the measurement names
    dist.measurements <- specimen.measurements[!grepl("\\.nm$", specimen.measurements)]
    angle.measurements <- specimen.measurements[grepl("\\.nm$", specimen.measurements)]
    
    # 3. Sample and assign (without replacement)
    perm.dist$long.name <- sample(dist.measurements, size = nrow(perm.dist))
    perm.angle$long.name <- sample(angle.measurements, size = nrow(perm.angle))
    
    # 4. Combine back into one dataframe
    perm.predictions <- bind_rows(perm.dist, perm.angle)
    
    ###
    
    result <- tryCatch({
      compare_specimen_percentiles(
        data = data,
        specimen = specimen,
        predictions = perm.predictions,
        maxNApercent = maxNApercent,
        normalize = normalize,
        calc.pval = calc.pval,
        norm_method = norm_method,
        remove.outliers = remove.outliers,
        percentile.dist = percentile.dist,
        cor.adj = cor.adj,
        adj.m = adj.m
      )
    }, error = function(e) {
      cat("Permutation failed:", conditionMessage(e), "\n")
      return(NULL)
    })
    
    if (is.null(result)) next
    
    perm_result <- result[[1]]
    perm_result$sampled_predictions <- paste(perm.predictions$long.name, collapse = ',')
    per.sum.abs.d <- c(per.sum.abs.d, perm_result$sum.abs.d)
    
    relevant.perm.count <- relevant.perm.count + 1
    permutations[relevant.perm.count, ] <- perm_result
  }
  
  return(permutations[1:relevant.perm.count, ])
}









## plot.group.as.specimen 4, now with Wilcoxon Rank Sum U test and added 'vs' parameter
plot.group.as.specimen4 = function(data, specimen,predictions,normalize = FALSE,
                                   modern = FALSE,norm_method = 'pensize',
                                   maxNApercent = 100,dir,calc.pval = 'U-test',
                                   mode = 'specimen',vs = 'both', output.dir,
                                   percentile.dist,remove.outliers = TRUE,
                                   cor.adj = TRUE,adj.m = c()){
  

  
  #### get.per Function ####
  # Get Percentile value of a specified specimen compared with a group (AMH,NE,ER)
  # Can normalize values if required
  # Required parameters: f.data, specific measurement,prediction data, pensize data
  get.per = function(m,specimen,group,modern = FALSE,remove.outliers = remove.outliers){
    
    #skip measurement if it does not exist in the specified specimen
    if(specimen %in% m$specimen == FALSE){
      cat('\t\tN/A\n')
      return(NA)
    }
    
    ## Get test value and compare values
    test.value = m[m$specimen == specimen,'value']
    cat(paste('\t\tvalue:',round(test.value,3)),'\n',sep='')
    
    #Remove test subject from the data
    m = m[m$specimen != specimen,]
    if(group == 'AMH' & modern == FALSE){
      compare.values = m[m$group %in% c('UPS','EHS'),'value']
    }else if(group == 'AMH' & modern == TRUE){
      compare.values = m[m$group %in% c('UPS','EHS','recentAMH'),'value']
    } else{
      compare.values = m[m$group == group,'value']
    }

    if (sum(!is.na(compare.values)) < 5) {
      cat("Skipping iteration due to insufficient non-NA values\n")
      return(NA)
    }

    if(remove.outliers ==TRUE){
      compare.values = compare.values[!compare.values %in% boxplot.stats(compare.values)$out] #removes outliers
    }
    
    if(percentile.dist == 'quantile'){
      
      compare.values = compare.values[order(compare.values)]
      N = length(compare.values)
      Quant = (1:N - 0.5)/N

      group.cdf = approxfun(compare.values,Quant,yleft = 0,yright = 1, ties = 'mean') #for interpolation
      
      percentile = 100*group.cdf(test.value)
      return(percentile)
    }
  }
  
  #### Gets all the relevant percentiles for the analysis ####
  cat(paste('\n',specimen,'\n'))
  groups = c('AMH','NE')
  # Initiate percentiles dataframe
  percentiles = data.frame(matrix(ncol = 6,nrow = 0))
  names(percentiles) = c('related.prediction','prediction','Vs','direction','percentile','loc')
  
  #get specimen group
  if(mode == 'specimen'){
    specimen.group = data[data$specimen == specimen,'group']
  }else if(mode == 'group'){
    specimen.group = specimen
  }

  relevant.measurements = predictions%>%
    filter(include.in.analysis == 'yes' & loc %in% c(1,2))%>%
    pull(long.name)
  
  
  #### Normalize if required ####
  
  mList = list()
  for(measurement in relevant.measurements){
    fixed.measurement = measurement
    cat(paste('\t',fixed.measurement,'\n'),sep='')
    measurement = gsub('[() :,/\\-]','.',measurement)
    m = get.m(data,fixed.measurement,predictions,normalize = normalize,
              method = norm_method,TS = specimen.group,adj.m = adj.m, cor.adj = cor.adj)
    mList[[measurement]] = m
    
    if(mode == 'group'){
      test.group = c(specimen,specimen,mean(m[m$group == specimen,'value']))
      m = rbind(m,test.group)
      m$value = as.numeric(m$value)
    }
    
    for(group in groups){
      
      # Skip if there aren't enough comparison values
      percentile <- tryCatch({
        get.per(m,specimen,group,modern = modern,remove.outliers = remove.outliers)
        }, error = function(e) {
        browser()
        percentile =   get.per(m,specimen,group,modern = modern,remove.outliers = remove.outliers)
        return(NA)  # or NA or some fallback
      })
      
      
      percentiles[nrow(percentiles)+1,] =
        c(related.prediction = as.character(predictions[predictions$long.name == fixed.measurement,'related.prediction']),
          measurement,
          group,
          direction = as.character(predictions[predictions$long.name == fixed.measurement,group]),
          percentile,
          loc =predictions[predictions$long.name == fixed.measurement,'loc'])
    }
  }
  
  percentiles[percentiles$prediction == 'Maximum.biparietal.breadth..Rightmire.et.al...2006...Cranial.vault..BP','prediction'] = 'Max biparietal braedth'

  #### Prepare percentiles for plot (discrete X to continous X for convenience) ####
  ## Initiate direction to graph location library
  dir.loc = hash()
  dir.loc[['less']] = 1
  dir.loc[['undetected']] = 2
  dir.loc[['greater']] = 3
  ## Replace direction with location
  percentiles$direction = values(dir.loc, keys = percentiles$direction)
  
  ## Remove NAs
  percentiles = percentiles[!is.na(percentiles$percentile),]
  percentiles$percentile = as.numeric(percentiles$percentile)
  
  ## Combine several measurements for the same prediction
  percentiles$prediction = NULL

  agg.perc = aggregate(percentile~.,percentiles,FUN=mean)
  agg.perc = subset(agg.perc,direction !=2)

  
  agg.perc = agg.perc[order(agg.perc$Vs,decreasing = TRUE),]
  # to remove AMH when AMH and NE have the exact same prediction
  
  
  if(vs == 'AMH'){
    agg.perc = subset(agg.perc,Vs == 'AMH')
  }else if(vs == 'NE'){
    agg.perc = subset(agg.perc,Vs == 'NE')
  }else if(vs == 'both.strict'){
    # to remove AMH when AMH and NE have the exact same prediction
    agg.perc = agg.perc[!duplicated(agg.perc[,c(1,3)]),]
    
  }

    ## Calculate total P.val
  if(!is.na(calc.pval)){
    agg.perc$percentile = agg.perc$percentile/100
    match = ''
    for(pred in 1:nrow(agg.perc)){
      if(agg.perc[pred,'direction'] == 3){
        agg.perc[pred,'d'] = agg.perc[pred,'percentile']-0.5
      }else if(agg.perc[pred,'direction'] == 1){
        agg.perc[pred,'d'] = (1-agg.perc[pred,'percentile']) - 0.5
      }
    }
    #Print all percentiles
    print.data.frame(agg.perc)
    
    
    n = nrow(agg.perc)
    mean.perc = round(mean(agg.perc$d),3)
    n.matched = sum(agg.perc$d>0)
    
    EKB = data[data$specimen ==specimen,'EKB.vis']
    if(calc.pval == 'mid.range.binomial'){
      x = sum(agg.perc$d>0)
      n = nrow(agg.perc)
      if(x==n){
        binom.test = binom.test(x,n,p = 0.5,alternative = 'g')
        binom.p.val=0.5*round(binom.test$p.value,7)
      }else{
        binom.test = binom.test(x,n,p = 0.5,alternative = 'g')
        binom.p.val1 = round(binom.test$p.value,7)
        
        binom.test = binom.test(x+1,n,p = 0.5,alternative = 'g')
        binom.p.val2 = round(binom.test$p.value,7)
        
        binom.p.val = 0.5*(binom.p.val1-binom.p.val2) +binom.p.val2
      }
      
      wilcox.test = wilcoxsign_test(agg.perc$d ~ rep(0,length(agg.perc$d)),
                                    distribution = "exact", alternative = 'g')
      wilcox.p.val = round(wilcox.test@distribution@pvalue(wilcox.test@statistic@teststatistic),7)
    }
    
    greater_than_AMH = sum(agg.perc$Vs == 'AMH'& agg.perc$direction == 3)
    less_than_AMH = sum(agg.perc$Vs == 'AMH'& agg.perc$direction == 1)
    greater_than_NE = sum(agg.perc$Vs == 'NE'& agg.perc$direction == 3)
    less_than_NE = sum(agg.perc$Vs == 'NE'& agg.perc$direction == 1)

    # Permutate measurements to coun
    sum.abs.d = sum(abs(agg.perc$d))

    pval = data.frame('specimen'=specimen,
                  'group' = as.character(specimen.group),
                  'binom.p.val' = binom.p.val,
                  'binom.score' = -log10(binom.p.val),
                  'wilcox.p.val' = wilcox.p.val,
                  'wilcox.score' = -log10(wilcox.p.val),
                  'comb_score' = NA,
                  'mean.d' = mean.perc,
                  'sum.abs.d' = sum.abs.d,
                  'total_predictions' = n,
                  'match_d_over_0.5' = n.matched,
                  'match_rate' = round(n.matched/n,3),
                  'perm_success_binom' = NA,
                  'perm_success_matched_d' = NA,
                  'perm_success_U' = NA,
                  'perm_success_mean_d' = NA,
                  'perm_success_dist' = NA,
                  'permutations'= NA,
                  'available.measurements'= NA,
                  'normalize'= normalize,
                  'norm_method' = norm_method,
                  'EKB' = EKB,
                  'recentAMH' = modern,
                  'greater_than_AMH' = greater_than_AMH,
                  'less_than_AMH' = less_than_AMH,
                  'greater_than_NE' = greater_than_NE,
                  'less_than_NE' = less_than_NE)
    
    pval$comb_score = sqrt(pval$binom.score^2+pval$wilcox.score^2)
    
    return(list(pval,mList,agg.perc))
  }
  
}


## get.m from 30-08-2022, takes relevant data and prepares it for plot
# needed for plot.prediction and plot.group
get.m = function(data,measurement,predictions,normalize = FALSE,TS = NA,
                 method = 'pensize', modern = FALSE,maxNApercent = 100,adj.m,cor.adj){

  ## Preprocess
  fixed.measurement = gsub('[(): ,/\\-]','.',measurement)

  prediction = predictions[predictions$long.name == fixed.measurement,]
  
  data = filter(data, group %in% c('UPS','EHS','NE','MPH','LMPEA','ERC_MP','ERC_EP'))
  
  m = select(data, c('specimen','group',measurement,'EKB.vis','Cranial.capacit'))%>%
    rename(value = 3)%>%
    filter(!is.na(value))
  
  #names(m)[3] = c('value')
  
  
  Neurocranium <- data %>%
    select(ends_with(".neu")) %>%
    colnames()
  
  Viscerocranium <- data %>%
    select(ends_with(".vis")) %>%
    colnames()
  
  cranial.loc = prediction$loc
  m.type = prediction$type
  
  if (m.type == 'distance') {
    if (normalize) {
      if (method == 'pensize') {
        norm.measurements = switch(
          cranial.loc, Neurocranium, Viscerocranium)  # , mandible, dentition
        
        norm.measurements = append(colnames(data)[1:7], norm.measurements)
        m.for.pensize = colnames(data)[colnames(data) %in% norm.measurements]
        m.for.pensize = m.for.pensize[m.for.pensize != measurement]
        norm.data = data[, colnames(data) %in% m.for.pensize]
        
        # Calculate and add pensize to m
        pensize = get.pensize3(norm.data, maxNApercent = maxNApercent) %>%
          select(specimen, pensize)
        
        m = m %>%
          left_join(pensize, by = 'specimen', copy = FALSE, keep = FALSE)
        
        # Standardize values
        m$value = scale(m$value) - m$pensize
        
      } else if (method == 'Cranial.capacit') {
        cat(paste('\t\t value was normalized by CC\n'))
        m$Cranial.capacit <- m$Cranial.capacit^(1/3)
        m$value = m$value / m$Cranial.capacit
        m = m[!is.na(m$value), ]
        
        
      } else if (method == 'biorbital_breadth') {
        m$value = m$value / m$EKB
        m = m[!is.na(m$value), ]
      }
    }
    
  } else {
    print('Measurement does not require normalization')
  }
  
  # Keep only key columns
  m = m[, 1:3]
  
  # Always apply correlation adjustment at the end
  if (cor.adj == TRUE & fixed.measurement %in% adj.m) {
    m$value = -m$value
  }
  
  return(m)
}
  

##Replace specimen names for graphs:
fix.names = function(names){
  
  og=c("irhoud_1",
       "irhoud_2",
       "irhoud_dentation_and_mandible",
       "florisbad",
       "omo_ii",
       "lh_18",
       "skhul_v",
       "skhul_ix",
       "qafzeh_ix",
       "qafzeh_5_6_7_11_15_27",
       "mladec_i",
       "mladec_ii",
       "mladec_v",
       "mladec_vi",
       "cro_magnon_i",
       "cro_magnon_ii",
       "cro_magnon_iii",
       "oase",
       "upper_cave_101",
       "upper_cave_103",
       "liujiang",
       "ate9_1",
       "atd6_15_atd6_69_atd6_96",
       "Antecessor_dental",
       "peking_x",
       "peking_xii",
       "peking_xiii",
       "peking_lii",
       "peking_rc_1996",
       "peking_dental",
       "nanjing1",
       "hexian",
       "sambungmacan_1",
       "sambungmacan_3",
       "sangiran_2",
       "sangiran17",
       "sangiran_dental",
       "ngandong_7",
       "ngandong_9",
       "ngandong_12",
       "dmanisi_211_2282",
       "dmanisi_2280",
       "dmanisi_2700_2735",
       "dmanisi_4500_2600",
       "sale",
       "knm_wt_15000",
       "er_3733",
       "er_3883",
       "oh24",
       "oh7",
       "er1805",
       "mauer_1",
       "arago_xxi_xlvii",
       "arago_xiii",
       "arago_ii",
       "broken_hill",
       "petralona_1",
       "ceprano",
       "steinheim_s11",
       "narmada",
       "eliye_springs",
       "ndutu",
       "Homo_sp._rabat",
       "Homo_sp._stw53",
       "oh_9",
       "sh4",
       "sh5",
       "tabun_c1",
       "tabun_c2",
       "spy_i",
       "spy_ii",
       "gibraltar1",
       "amud",
       "la_chapelle_aux_saint",
       "la_ferrassie_1",
       "shanidar_1",
       "shanidar_2",
       "shanidar_5",
       "cesaire",
       "saccopastore_i",
       "saccopastore_ii",
       "Neanderthal_type",
       "saldanha",
       "bodo",
       "ternifine_1",
       "ternifine_2",
       "ternifine_3",
       "ternifine_4",
       "xiahe",
       "Homo_altaiensis",
       "dali",
       "hualongdong",
       "Homo_longi",
       "jinniushan",
       "maba",
       "xuchang")
  
  new=c("Jebel Irhoud 1",
        "Jebel Irhoud 2",
        "Jebel Irhoud Dentition and Mandible",
        "Florisbad",
        "Omo II",
        "LH 18",
        "Skhul V",
        "Skhul IX",
        "Qafzeh IX",
        "Qafzeh 5 6 7 11 15 27",
        "Mladec I",
        "Mladec II",
        "Mladec V",
        "Mladec VI",
        "Cro-Magnon I",
        "Cro-Magnon II",
        "Cro-Magnon III",
        "Oase",
        "Zhoukouodian UP 101",
        "Zhoukouodian UP 103",
        "Liujiang",
        "ate9 1",
        "atd6 15 atd6 69 atd6 96",
        "Antecessor Dental",
        "Peking X",
        "Peking XII",
        "Peking XIII",
        "Peking LII",
        "Peking Reconstrction",
        "Peking Dental",
        "Nanjing 1",
        "Hexian",
        "Sambungmacan 1",
        "Sambungmacan 3",
        "Sangiran 2",
        "Sangiran 17",
        "Sangiran Dental",
        "Ngandong 7",
        "Ngandong 9",
        "Ngandong 12",
        "Dmanisi 211 2282",
        "Dmanisi 2280",
        "Dmanisi 2700 2735",
        "Dmanisi 4500 2600",
        "Sale",
        "KNM-WT 15000",
        "KNM-ER 3733", #Kenya National Museum -East Rudolf 
        "KNM-ER 3883",
        "OH 24",
        "OH 7",
        "KNM ER 1805",
        "Mauer 1",
        "Arago XXI XLVII",
        "Arago XIII",
        "Arago II",
        "Kabwe 1",
        "Petralona 1",
        "Ceprano",
        "Steinheim",
        "Narmada",
        "Eliye Springs",
        "Ndutu",
        "Rabat",
        "STW 53",
        "OH 9",
        "Sima de Los Huesos 4",
        "Sima de Los Huesos 5",
        "Tabun C1",
        "Tabun C2",
        "Spy I",
        "Spy II",
        "Gibraltar 1",
        "Amud",
        "La Chapelle aux Saints 1",
        "La Ferrassie 1",
        "Shanidar 1",
        "Shanidar 2",
        "Shanidar 5",
        "Cesaire",
        "Saccopastore I",
        "Saccopastore II",
        "Neanderthal 1",
        "Saldanha",
        "Bodo",
        "Ternifine 1",
        "Ternifine 2",
        "Ternifine 3",
        "Ternifine 4",
        "Xiahe",
        "Denisovan",
        "Dali",
        "Hualongdong",
        "Harbin",
        "Jinniushan",
        "Maba",
        "Xuchang")
  
  new_vec = stringr::str_replace_all(names,setNames(new,og))
  
  
  return(new_vec)
}


