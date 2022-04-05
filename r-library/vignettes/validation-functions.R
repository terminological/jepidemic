#options(warn=2)
#options(warn=1)

## Lag analysis ----

#' This does a lag analysis on the sawtooth timeseries with a fixed infectivity profile.
#'
#' @param estimators - a tibble of estimator methods to compare. This should have 2 columns - "model" the model name, and estimFn: the estimation function as a list column.
#' This function should take 2 parameters, "ts" which will be a simple time series and "infectivityProfile" which will be the passed the infectivity profile yMatrix  
#' @param trueVar - the column in the ts time series which contains the actual value of the growth rate or reproduction number
#' @param estimateVar - the column in the estimated time series that we can compare to the true value of the the estimate.
#'
#' @return a list containing 3 elements: the "estimates" containing the estimates from each of the methods being compared, the "rmseByOffset" which summarises the
#' estimates for each model compared to the actual value with a series of offsets, and the "modelLag", the value of the offset woth the minimum rmse of each of the methods
#' being compared.
#' @export 
lagAnalysis = function(estimators, estimation = "Rt") {
  trueVar = as.symbol(paste0(estimation, ".actual"))
  estimateVar = as.symbol(paste0(estimation, ".Quantile.0.5"))
  
  triangularSim = jepidemic::lagAnalysisDataset
  #triangular = enframe(triangularSim) %>% pivot_wider()
  triangular = tibble(name = triangularSim$name, ts = list(triangularSim$ts), infectivityProfile = list(triangularSim$infectivityProfile))
  
  periodicEstimates = triangular %>% inner_join(estimators, by=character()) %>% mutate(
    
    estimate = pmap(.l=list(fn = estimFn, ts = ts, infectivityProfile = infectivityProfile), .f = function(fn,ts,infectivityProfile) {
      fn(ts,infectivityProfile)
    })
  )

  periodicRt = periodicEstimates %>% select(model,estimate) %>% unnest(estimate)
  periodicTs = periodicEstimates %>% select(model,ts) %>% unnest(ts)

  offset = tibble(shift = c(-7:21)) %>% inner_join(
    periodicTs %>% select(theoretical.date = date, !!trueVar), by=character()) %>% 
    mutate(
      estimate.date = theoretical.date+shift
    )
  
  offsetPeriodicRt = periodicRt %>% select(-any_of(as_label(trueVar))) %>% 
    inner_join(offset, by=c("date"="estimate.date"), suffix=c("",".offset")) %>%
    filter(date >= min(date)+28 & date <= max(date)-28)
  
  offsetSummary = offsetPeriodicRt %>% 
    mutate(sqrdErr = (!!trueVar-!!estimateVar)^2) %>%
    group_by(model,shift) %>%
    summarise(rmse = sqrt(mean(sqrdErr,na.rm = TRUE)))
  
  offsetSummarySummary = offsetSummary %>% group_by(model) %>% group_modify(function(d,g,...) {
    out = tryCatch({
      lagFn = splinefun(d$shift, d$rmse)
      min = d$shift[d$rmse==min(d$rmse,na.rm = TRUE)]
      root = uniroot(lagFn, lower = min-2, upper = min+2, deriv=1)
      return(tibble(medianLag = root$root))
    }, error = function(e) {
      return(tibble(medianLag = d$shift[d$rmse==min(d$rmse,na.rm = TRUE)]))
    })
    return(out)
  }) %>% mutate(label = sprintf("%1.2f days",medianLag))

  out = list(
    estimates = periodicRt %>% inner_join(periodicTs, by=c("model","date"), suffix=c(".estimate","")),
    rmseByOffset = offsetSummary,
    modelLag = offsetSummarySummary
  )
  
  class(out) = "lag_analysis_result"
  
  return(out)
}

## Lag analysis visualisation ----

lagTimeseriesPlot = function(lagAnalysisResult, estimation = "Rt", ylim=c(0.4,1.6), ncol=1) {
  estimates = lagAnalysisResult$estimates
  trueVar = as.symbol(paste0(estimation, ".actual"))
  estimateVar = as.symbol(paste0(estimation, ".Quantile.0.5"))
  lowVar = as.symbol(paste0(estimation, ".Quantile.0.025"))
  highVar = as.symbol(paste0(estimation, ".Quantile.0.975"))
  
  suppressWarnings(
    ggplot(estimates %>% filter(!is.nan(!!estimateVar)), mapping = aes(x=date))+
      geom_ribbon(aes(ymin = !!lowVar, ymax = !!highVar), alpha=0.2,colour = NA,fill="grey20")+
      geom_line(aes(y=!!estimateVar),colour="black")+
      guides(colour="none")+
      geom_line(mapping = aes(x=date,y=!!trueVar), inherit.aes = FALSE, colour="red")+
      xlab(NULL)+
      ylab(paste0(estimation," estimate"))+
      facet_wrap(vars(model),ncol=ncol)+coord_cartesian(ylim=ylim)
  )
}

lagOffsetPlot = function(lagAnalysisResult, ncol=1) {
  labelSize = 8/ggplot2:::.pt/(96/72)
  ggplot(lagAnalysisResult$rmseByOffset, aes(x=shift,y=rmse))+
    geom_point()+
    geom_line()+
    geom_vline(data = lagAnalysisResult$modelLag, mapping=aes(xintercept = medianLag), colour="blue")+
    geom_label(data = lagAnalysisResult$modelLag, mapping=aes(label=label,x = medianLag), y=Inf, fill=NA, label.size=0, label.padding = unit(1, "lines"),hjust=0, vjust="inward", colour="blue",size=labelSize)+
    ylab("RMSE")+xlab("shift (days)")+
    facet_wrap(vars(model),ncol=ncol)
}

lagPlot = function(lagAnalysisResult, ...) {
  p1 = lagTimeseriesPlot(lagAnalysisResult, ...)
  p2 = lagOffsetPlot(lagAnalysisResult, ncol=1)
  p3 = p1+p2+plot_layout(nrow=1,widths = c(2,1))+plot_annotation(tag_levels = "A")
  return(p3)
}


# reconstruct a CDF on a regular grid from the quantiles
inferQuantiles = function(rtEstimate, prefix="Rt") {
  nm = function(s) as.symbol(paste0(prefix,s))
  quantilePrefix = paste0(prefix,".Quantile.")
  rtEstimate %>% 
    mutate(
      median = !!nm(".Quantile.0.5"),
      lo = !!nm(".Quantile.0.025"),
      hi = !!nm(".Quantile.0.975"),
      iqr = !!nm(".Quantile.0.75")-!!nm(".Quantile.0.25")
    ) %>%
    pivot_longer(cols = starts_with(quantilePrefix), names_to = "q", values_to = "x") %>% 
    mutate(q = as.numeric(stringr::str_remove(q,quantilePrefix))) %>% 
    filter(!is.na(x)) %>%
    nest(quantile.data = c(q,x)) %>%
    mutate(quantile.model = map(quantile.data, function(data) {
      range = max(data$x,na.rm = TRUE)-min(data$x,na.rm = TRUE)
      extendedData = c(
        -100,
        min(data$x,na.rm = TRUE)-range/2,
        data$x,
        max(data$x,na.rm = TRUE)+range/2,
        100)
      splinefun(extendedData, c(0,0,data$q,1,1), method = "monoH.FC")
    })) %>%
    mutate(quantile.model.range = map(quantile.data, function(data) {
      range = max(data$x,na.rm = TRUE)-min(data$x,na.rm = TRUE)
      extendedData = c(
        min(data$x,na.rm = TRUE)-range/2,
        max(data$x,na.rm = TRUE)+range/2
      )
      return(extendedData)
    })) #%>%
  # mutate(inv.quantile.model = map(quantile.data, function(data) {
  #   range = max(data$x,na.rm = TRUE)-min(data$x,na.rm = TRUE)
  #   extendedData = c(
  #     min(data$x,na.rm = TRUE)-range/2,
  #     data$x,
  #     max(data$x,na.rm = TRUE)+range/2)
  #   splinefun(c(0,data$q,1), extendedData, method = "monoH.FC")
  # })) 
}

calculateMetrics = function(quantileEstimate, originalTs, medianLag, estimation = "Rt", criticalValue = 1) {
  
  trueVar = as.symbol(paste0(estimation, ".actual"))
  estimateVar = as.symbol(paste0(estimation, ".Quantile.0.5"))
  lowVar = as.symbol(paste0(estimation, ".Quantile.0.025"))
  highVar = as.symbol(paste0(estimation, ".Quantile.0.975"))
  
  trueVar = ensym(trueVar)
  originalTs %>% 
    select(bootstrap,date, actual = !!trueVar) %>%
    left_join( # here is the point where we truncate the data set as a result of shifting the estiamtes by the lag. If this were a left join
      quantileEstimate %>% 
        mutate(
          estimate.date = date,
          date = estimate.date-round(medianLag)),
      by=c("bootstrap","date")
    ) %>%
    mutate(
      boundaryEffect = case_when(
        date < min(date)+21 ~ "start",
        date >= max(date)-21 ~ "end",
        TRUE ~ "middle"
      ) %>% ordered(levels=c("start","middle","end"))
    ) %>%
    mutate(
      quantile.actual = map2_dbl(quantile.model, actual, function(.x,.y) if(is.null(.x)) return(NA_real_) else return(.x(.y))),
      residual = median-actual,
      calibration = ifelse(actual > lo & actual < hi,1,0),
      critical_threshold = ifelse(sign(actual-criticalValue)!=sign(lo-criticalValue) & sign(actual-criticalValue)!=sign(hi-criticalValue),1,0)
    ) %>%
    mutate(
      quantile.actual = case_when(
        quantile.actual < 0 ~ 0,
        quantile.actual > 1 ~ 1,
        TRUE ~ quantile.actual
      )
    ) %>% 
    mutate(
      crps = pmap_dbl(
        .l = list(F = quantile.model, limits = quantile.model.range, y = actual), 
        .f = function(F,limits,y) {
          if(is.null(F)) return(NA_real_)
          tmp = integrate(f = function(x) (F(x)-(x >= y))^2, lower = limits[1], upper=limits[2])
          return(tmp$value)
        }
      )
    )
}

qualityAnalysis = function(estimators, modelLag, estimation = "Rt", criticalValue = 1, synthetic = jepidemic::validationDataset) {
  
  syntheticEstimates(estimators,synthetic) %>%
    validationMetrics(modelLag, estimation)
  
}

syntheticEstimates = function(estimators, synthetic = jepidemic::validationDataset) {
  syntheticEstimates = synthetic %>% inner_join(estimators, by=character()) %>% 
    mutate(
      estimate = pmap(.l=list(fn = estimFn, ts = ts, infectivityProfile = infectivityProfile), .f = function(fn,ts,infectivityProfile) {
        ts %>% group_by(bootstrap) %>% group_modify(function(d,g,...) {
          fn(d, infectivityProfile)
        })
      })
    )
}

validationMetrics = function(syntheticEstimates,modelLag, estimation = "Rt", criticalValue = 1) {
  syntheticEstimates %>% 
    mutate(
      quantileEstimate = map(estimate, inferQuantiles, prefix=estimation)
    ) %>%
    inner_join(modelLag, by="model") %>%
    mutate(
      comparison = pmap(.l=list(quantileEstimate, ts, medianLag), .f=calculateMetrics, estimation=estimation, criticalValue = criticalValue)
    ) %>%
    select(-quantileEstimate)
} 

## Divergence and comparison metrics ----
# KL divergence for 
# 
# klGammaVsGamma = function(mean1, sd1, mean2, sd2) {
#   if (sd1>mean1 | sd2>mean2) stop("sd must be smaller than mean")
#   range = seq(0,10,length=200)
#   g1 = dgamma(range, shape=mean1^2/sd1^2,rate = mean1/sd1^2)
#   g2 = dgamma(range, shape=mean2^2/sd2^2,rate = mean2/sd2^2)
#   kl = flexmix::KLdiv(cbind(g1=g1,g2=g2))[1,2]
#   # A matrix of KL divergences where the rows correspond to using the respective distribution as  in the formula above.
#   return(kl)
# }
# 
# klGammaVsUnif = function(value1, precision1, mean2, sd2) {
#   if (sd2>mean2) stop("sd must be smaller than mean")
#   range = seq(0,10,length=200)
#   g1 = dunif(range, min = value1-precision1,max = value1+precision1)
#   g2 = dgamma(range, shape=mean2^2/sd2^2,rate = mean2/sd2^2)
#   kl = flexmix::KLdiv(cbind(g1=g1,g2=g2))[1,2]
#   # A matrix of KL divergences where the rows correspond to using the respective distribution as  in the formula above.
#   return(kl)
# }
# 
# 
# klGammaVsGamma(mean1=2,sd1=1.5,mean2=2.3,sd2=2.1)
# klGammaVsUnif(value1=2,precision1=0.0001,mean2=2,sd2=0.002)

## Validation plots ----

estimationExamplePlot = function(qualityAnalysisResult, estimation = "Rt", ylim = c(0.6,1.4), simFilterExpr = {weekendEffect == 0.1 & seed == 100}, modelFilterExpr = TRUE, bootstraps = 1) {
  
  trueVar = as.symbol(paste0(estimation, ".actual"))
  estimateVar = as.symbol(paste0(estimation, ".Quantile.0.5"))
  lowVar = as.symbol(paste0(estimation, ".Quantile.0.025"))
  highVar = as.symbol(paste0(estimation, ".Quantile.0.975"))
  
  simFilterExpr = enexpr(simFilterExpr)
  modelFilterExpr = enexpr(modelFilterExpr)
  
  tmp = qualityAnalysisResult %>% 
    filter(!!simFilterExpr) %>%
    filter(!!modelFilterExpr) %>%
    mutate(label = paste0(smoothLabel,"; ",seedLabel,"; ",weekendLabel))
  
  tmp2 = tmp %>% ungroup() %>% select(label,ts) %>% distinct() %>% unnest(ts) %>% filter(bootstrap %in% bootstraps) %>% mutate(row="simulation")
  tmp3 = tmp %>% ungroup() %>% select(model,label,estimate) %>% unnest(estimate) %>% filter(bootstrap %in% bootstraps)
  
  p1 = ggplot(tmp2, aes(x=date,y=value,group=bootstrap))+
    geom_point(size=0.5)+
    geom_line(alpha=0.1)+
    ylab("cases")+
    geom_line(aes(y=Est.actual),colour="red")+facet_grid(rows = vars(row), cols=vars(label))
  
  p2 = ggplot(tmp3 %>% filter(!is.nan(!!estimateVar)), mapping = aes(x=date, group=bootstrap))+
    geom_ribbon(aes(ymin = !!lowVar, ymax = !!highVar), alpha=0.2,colour = NA,fill="grey20")+
    geom_line(aes(y=!!estimateVar)) +
    guides(colour="none")+
    coord_cartesian(ylim=ylim)+
    geom_line(data=tmp2,mapping=aes(x=date,y=!!trueVar),colour="red")+
    ylab(paste0(estimation," estimate"))+xlab(NULL) +facet_grid(rows = vars(model),cols = vars(label))
  
  p3 = p1+standardPrintOutput::hideX() + p2 +patchwork::plot_layout(ncol=1, heights=c(1,4))+patchwork::plot_annotation(tag_levels = "A")
  return(p3)
}


percentage = function(groupedDf,var) {
  var = ensym(var)
  groupedDf %>% summarise(
      tmp_sum = sum(!!var,na.rm=TRUE), N = sum(!is.na(!!var)),
      .groups="drop"
    ) %>%
    mutate(middle = tmp_sum/N)
}

boxplotData = function(groupedDf, var) {
  var = ensym(var)
  groupedDf %>% summarise(
      tibble(
        names = c("ymin","lower","middle","upper","ymax"),
        values = boxplot.stats(!!var, do.conf = FALSE, do.out = FALSE)$stats,
      ), .groups="drop_last"
    ) %>% pivot_wider(names_from = names, values_from = values) %>%
    ungroup()
}

densityData = function(groupedDf, var,bw=0.025) {
  var = ensym(var)
  groupedDf %>% summarise(
    tibble(
      x=seq(0,1,length.out = 1001),
      y=density(x = !!var, from = 0, to = 1, n = 1001, na.rm = TRUE,bw=bw)$y,
    ) %>% mutate(
      dydx = signal::sgolayfilt(y, n = 11, p=2,ts = 1/1000,m = 1)
    ), 
    .groups="drop_last"
  )
}

estimateSummaryPlot = function(qualityAnalysisResult, errorLimits = NULL, crpsLimits = NULL, critLimits = NULL) {
  tmp6 = qualityAnalysisResult %>% 
    ungroup() %>%
    select(model,weekendLabel,seedLabel,smoothLabel,comparison) %>% 
    unnest(comparison) 
  
  d1 = tmp6 %>% group_by(model) %>% boxplotData(residual)
  d2 = tmp6 %>% group_by(model) %>% percentage(calibration)
  d3 = tmp6 %>% filter(!is.na(quantile.actual))
  d4 = tmp6 %>% group_by(model) %>% boxplotData(crps)
  d5 = tmp6 %>% group_by(model) %>% percentage(critical_threshold)
  
  lim = max(abs(c(
    d1 %>% pull(ymax),
    d1 %>% pull(ymin)
  )))
  if(is.null(errorLimits)) errorLimits = c(-lim*1.05,lim*1.05)
  calibrationLimits = c(0,1)
  if(is.null(crpsLimits)) crpsLimits = c(0,max(d4 %>% pull(ymax),na.rm = TRUE))
  if(is.null(critLimits)) critLimits = c(0,max(d5 %>% pull(middle),na.rm = TRUE))
  
  # residual error
  p1 = ggplot(d1, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=model, colour=model))+geom_boxplot(stat="identity",size=0.2)+coord_cartesian(ylim = errorLimits)+xlab(NULL)+ylab("bias")
  # calibration
  p2 = ggplot(d2,aes(x=model, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+coord_cartesian(ylim = calibrationLimits)+ylab("calibration")+xlab(NULL)
  # Rank histogram / Quantile density
  p3 = ggplot(d3, aes(y=quantile.actual, x=model, colour=model))+geom_violin(bw=0.025,size=0.2)+xlab(NULL)+ylab("quantile density")
  # Continuous ranked probability score
  p4 = ggplot(d4, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=model, colour=model))+geom_boxplot(stat="identity",size=0.2)+xlab(NULL)+ylab("CRPS")+coord_cartesian(ylim = crpsLimits)
  # Critical threshold
  p5 = ggplot(d5, aes(x=model, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+ylab("crit threshold")+xlab(NULL)+coord_cartesian(ylim = critLimits)
  
  p = p1+guides(colour=guide_none())+standardPrintOutput::hideX()+
    p2+standardPrintOutput::smallLegend(spaceLegend = 1)+guides(colour=guide_legend(title=NULL))+standardPrintOutput::hideX()+
    p3+guides(colour=guide_none())+standardPrintOutput::hideX()+
    p4+guides(colour=guide_none())+standardPrintOutput::hideX()+
    p5+guides(colour=guide_none())+standardPrintOutput::hideX()+
    patchwork::guide_area()+
    patchwork::plot_annotation(tag_levels = "A")+
    patchwork::plot_layout(ncol=3,guides = "collect")
  return(p)
}

estimateBreakdownPlot = function(qualityAnalysisResult, errorLimits = NULL, crpsLimits = NULL, critLimits = NULL) {
  tmp5 = qualityAnalysisResult %>% 
    ungroup() %>%
    select(model,weekendLabel,seedLabel,smoothLabel,comparison) %>% 
    unnest(comparison) 
  
  # residual error
  d1 = tmp5 %>% group_by(model,smoothLabel) %>% boxplotData(residual)
  d2 = tmp5 %>% group_by(model,weekendLabel) %>% boxplotData(residual)
  d3 = tmp5 %>% group_by(model,seedLabel) %>% boxplotData(residual)
  d4 = tmp5 %>% group_by(model,boundaryEffect) %>% boxplotData(residual)
  errorMax = bind_rows(d1,d2,d3,d4) %>% pull(ymax) %>% max(na.rm = TRUE)
  errorMin = bind_rows(d1,d2,d3,d4) %>% pull(ymin) %>% min(na.rm = TRUE)
  lim = max(abs(c(errorMax,errorMin)))
  if (is.null(errorLimits)) errorLimits = c(-lim*1.05,lim*1.05)
  
  p1 = ggplot(d1, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=smoothLabel, colour=model))+geom_boxplot(stat="identity",size=0.2)+coord_cartesian(ylim = errorLimits)+xlab(NULL)+ylab("error / bias")
  p2 = ggplot(d2, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=weekendLabel, colour=model))+geom_boxplot(stat="identity",size=0.2)+coord_cartesian(ylim = errorLimits)+xlab(NULL)+ylab("error / bias")
  p3 = ggplot(d3, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=seedLabel, colour=model))+geom_boxplot(stat="identity",size=0.2)+coord_cartesian(ylim = errorLimits)+xlab(NULL)+ylab("error / bias")
  p4 = ggplot(d4, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=boundaryEffect, colour=model))+geom_boxplot(stat="identity",size=0.2)+coord_cartesian(ylim = errorLimits)+xlab(NULL)+ylab("error / bias")
  
  # calibration
  calibrationLimits = c(0,1)
  p5 = ggplot(tmp5 %>% group_by(model,smoothLabel) %>% percentage(calibration),aes(x=smoothLabel, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+coord_cartesian(ylim = calibrationLimits)+ylab("calibration")+xlab(NULL)
  p6 = ggplot(tmp5 %>% group_by(model,weekendLabel) %>% percentage(calibration),aes(x=weekendLabel, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+coord_cartesian(ylim = calibrationLimits)+ylab("calibration")+xlab(NULL)
  p7 = ggplot(tmp5 %>% group_by(model,seedLabel) %>% percentage(calibration),aes(x=seedLabel, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+coord_cartesian(ylim = calibrationLimits)+ylab("calibration")+xlab(NULL)
  p8 = ggplot(tmp5 %>% group_by(model,boundaryEffect) %>% percentage(calibration),aes(x=boundaryEffect, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+coord_cartesian(ylim = calibrationLimits)+ylab("calibration")+xlab(NULL)
  
  # Rank histogram / Quantile density
  p9 = ggplot(tmp5 %>% filter(!is.na(quantile.actual)), aes(y=quantile.actual, x=smoothLabel, colour=model))+geom_violin(bw=0.025,size=0.2)+xlab(NULL)+ylab("quant. dens.")
  p10 = ggplot(tmp5 %>% filter(!is.na(quantile.actual)), aes(y=quantile.actual, x=weekendLabel, colour=model))+geom_violin(bw=0.025,size=0.2)+xlab(NULL)+ylab("quant. dens.")
  p11 = ggplot(tmp5 %>% filter(!is.na(quantile.actual)), aes(y=quantile.actual, x=seedLabel, colour=model))+geom_violin(bw=0.025,size=0.2)+xlab(NULL)+ylab("quant. dens.")
  p12 = ggplot(tmp5 %>% filter(!is.na(quantile.actual)), aes(y=quantile.actual, x=boundaryEffect, colour=model))+geom_violin(bw=0.025,size=0.2)+xlab(NULL)+ylab("quant. dens.")
  
  d13 = tmp5 %>% group_by(model,smoothLabel) %>% boxplotData(crps)
  d14 = tmp5 %>% group_by(model,weekendLabel) %>% boxplotData(crps)
  d15 = tmp5 %>% group_by(model,seedLabel) %>% boxplotData(crps)
  d16 = tmp5 %>% group_by(model,boundaryEffect) %>% boxplotData(crps)
  crpsMax = bind_rows(d13,d14,d15,d16) %>% pull(ymax) %>% max(na.rm = TRUE)
  if (is.null(crpsLimits)) crpsLimits = c(0,crpsMax*1.1)
  
  # Continuous ranked probability score
  p13 = ggplot(d13, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=smoothLabel, colour=model))+geom_boxplot(stat="identity",size=0.2)+xlab(NULL)+ylab("CRPS")+coord_cartesian(ylim = crpsLimits)
  p14 = ggplot(d14, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=weekendLabel, colour=model))+geom_boxplot(stat="identity",size=0.2)+xlab(NULL)+ylab("CRPS")+coord_cartesian(ylim = crpsLimits)
  p15 = ggplot(d15, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=seedLabel, colour=model))+geom_boxplot(stat="identity",size=0.2)+xlab(NULL)+ylab("CRPS")+coord_cartesian(ylim = crpsLimits)
  p16 = ggplot(d16, aes(ymin=ymin,lower=lower,middle=middle,upper=upper,ymax=ymax,x=boundaryEffect, colour=model))+geom_boxplot(stat="identity",size=0.2)+xlab(NULL)+ylab("CRPS")+coord_cartesian(ylim = crpsLimits)
  
  d17 = tmp5 %>% group_by(model,smoothLabel) %>% percentage(critical_threshold)
  d18 = tmp5 %>% group_by(model,weekendLabel) %>% percentage(critical_threshold)
  d19 = tmp5 %>% group_by(model,seedLabel) %>% percentage(critical_threshold)
  d20 = tmp5 %>% group_by(model,boundaryEffect) %>% percentage(critical_threshold)
  critMax = bind_rows(d17,d18,d19,d20) %>% pull(middle) %>% max(na.rm = TRUE)
  if (is.null(critLimits)) critLimits = c(0,critMax*1.1)
  
  p17 = ggplot(d17,aes(x=smoothLabel, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+ylab("crit. thresh.")+xlab(NULL)+coord_cartesian(ylim = critLimits)
  p18 = ggplot(d18,aes(x=weekendLabel, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+ylab("crit. thresh.")+xlab(NULL)+coord_cartesian(ylim = critLimits)
  p19 = ggplot(d19,aes(x=seedLabel, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+ylab("crit. thresh.")+xlab(NULL)+coord_cartesian(ylim = critLimits)
  p20 = ggplot(d20,aes(x=boundaryEffect, y=middle,colour=model))+geom_bar(stat="identity",fill="white",size=0.2,position="dodge")+ylab("crit. thresh.")+xlab(NULL)+coord_cartesian(ylim = critLimits)
  
  p = p1+standardPrintOutput::hideX()+guides(colour=guide_none())+
    p2+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p3+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p4+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p5+standardPrintOutput::hideX()+standardPrintOutput::smallLegend(spaceLegend = 0.75)+guides(colour=guide_legend(title=NULL))+
    p6+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p7+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p8+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p9+standardPrintOutput::hideX()+guides(colour=guide_none())+
    p10+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p11+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p12+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p13+standardPrintOutput::hideX()+guides(colour=guide_none())+
    p14+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p15+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p16+standardPrintOutput::hideX()+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p17+guides(colour=guide_none())+
    p18+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p19+standardPrintOutput::hideY()+guides(colour=guide_none())+
    p20+standardPrintOutput::hideY()+guides(colour=guide_none())+
    patchwork::plot_annotation(tag_levels = "A")+
    patchwork::plot_layout(ncol=4,guides = "collect")
  
  return(p)
}

estimateTimeseriesPlot = function(qualityAnalysisResult, simFilterExpr, modelFilterExpr, labelExpr = {paste0(smoothLabel,"; ",seedLabel,"; ",weekendLabel)}, errorLimits = NULL, crpsLimits = NULL) {
  simFilterExpr = enexpr(simFilterExpr)
  modelFilterExpr = enexpr(modelFilterExpr)
  labelExpr = enexpr(labelExpr)
  # This summary is on a per model basis
  # Here we filter out the specific model
  tmp = qualityAnalysisResult %>% 
    filter(!!modelFilterExpr) %>%
    filter(!!simFilterExpr) %>% 
    mutate(label = !!labelExpr)
  
  tmp4 = tmp %>% ungroup() %>% select(label,comparison) %>% unnest(comparison)
  tmp3 = tmp %>% ungroup() %>% select(label,ts) %>% unnest(ts)
  # spread skill analysis at level of individual estimate
  # thsi is is variance versus square of residual essentially.
  # ggplot(tmp4,aes(x=Rt.SD^2,y=Rt.residual.sqrd, colour=lambda_t))+geom_point()+facet_wrap(vars(label))
  
  # Rank histogram. Rank of actual versus bootstrap. Here as we have quantile we instead look at estimated quantile of actual.
  # Histogram shoudl be flat if confidence appropriate. If U shaped then too confident, if hill shaped then not confident enough.
  # ggplot(tmp4,aes(x=quantile.actual))+geom_density(bw=0.025)+facet_wrap(vars(label))
  
  # The theoretical values:
  p1 = ggplot(tmp3 %>% filter(bootstrap==1), aes(x=date))+
    geom_point(aes(x=date,y=value),size=0.1)+ 
    geom_line(aes(y=Est.actual), colour="red")+
    facet_wrap(vars(label),nrow = 1)+ylab(latex2exp::TeX("model cases"))+scale_y_continuous(trans="log1p",breaks = ukcovidtools::breaks_log1p())
  
  p2 = ggplot(tmp4 %>% filter(bootstrap==1), aes(x=date))+
    geom_ribbon(aes(ymin=lo, ymax=hi),fill="grey80",colour=NA)+
    geom_line(aes(y=median),colour="black")+
    geom_line(aes(y=actual), colour="red")+
    facet_wrap(vars(label),nrow = 1)+ylab(latex2exp::TeX("model $R_t$"))+
    coord_cartesian(ylim=c(0.6,1.4))
  
  # ggplot(tmp4, aes(x=date,y=quantile.residual))+geom_point(size=0.5)+plotRollingQuantiles(stepCDF3 %>% ungroup(), quantile.residual,window = 28)+facet_wrap(vars(label))
  lim = quantile(abs(tmp4$residual),probs = 0.975,na.rm = TRUE)
  if(is.null(errorLimits)) errorLimits = c(-lim,lim)
  
  # Timeseries of the comparison stats
  p3 = ggplot(tmp4 %>% filter(!is.na(residual)), aes(x=date,y=residual))+
    geom_point(size=0.1,alpha=0.25)+
    plotRollingQuantiles(tmp4 %>% group_by(label), residual,window = 14)+
    coord_cartesian(ylim=errorLimits)+ 
    ylab("bias")+
    facet_wrap(vars(label),nrow = 1)
  
  p4 = ggplot(tmp4, aes(x=date,y=calibration))+
    plotRollingProportion(tmp4 %>% group_by(label),calibration,window = 7) + 
    ylab("calibration")+
    facet_wrap(vars(label),nrow = 1)
  
  p5 = ggplot(tmp4, aes(x=date,y=quantile.actual))+
    plotRollingDensity(tmp4 %>% group_by(label), quantile.actual, window = 14) + 
    plotRollingDeciles(tmp4 %>% group_by(label), quantile.actual, window = 28) + 
    ylab("quantile density")+
    facet_wrap(vars(label),nrow = 1)+
    scale_fill_gradient(low="white",high="black",guide = "none")
  
  if(is.null(crpsLimits)) crpsLimits = c(0,quantile(abs(tmp4$crps),probs = 0.975,na.rm = TRUE))
  p6 = ggplot(tmp4, aes(x=date,y=crps))+
    geom_point(size=0.1,alpha=0.25)+
    plotRollingQuantiles(tmp4 %>% group_by(label), crps, window = 14) + 
    ylab("instantaneous CRPS")+facet_wrap(vars(label),nrow = 1)+
    coord_cartesian(ylim=crpsLimits)
  
  p7 = ggplot(tmp4, aes(x=date,y=critical_threshold))+
    plotRollingProportion(tmp4 %>% group_by(label), critical_threshold,window = 7) + 
    ylab("crit threshold")+facet_wrap(vars(label),nrow = 1)
  
  p = 
    p1+standardPrintOutput::hideX()+
    p2+standardPrintOutput::hideX()+
    p3+standardPrintOutput::hideX()+
    p4+standardPrintOutput::hideX()+
    p5+standardPrintOutput::hideX()+
    p6+standardPrintOutput::hideX()+
    p7+
    plot_annotation(tag_levels = "A")+plot_layout(ncol=1, heights = c(1,1,2,2,2,2,2))
  
  return(p)
}

printMedianAndCI = function(x) do.call(sprintf, c("%1.2g [95%% CI %1.2g \u2013 %1.2g]",as.list(quantile(x,c(0.5,0.025,0.975),na.rm = TRUE))))
printMedianAndCI3 = function(x) do.call(sprintf, c("%1.3f [95%% CI %1.3f \u2013 %1.3f]",as.list(quantile(x,c(0.5,0.025,0.975),na.rm = TRUE))))
printBinomCI = function(x) do.call(sprintf, c("%1.1f%% [95%% CI %1.1f%% \u2013 %1.1f%%]", as.list(binom::binom.confint(x=sum(x==1,na.rm=TRUE),n=length(na.omit(x)),methods = "wilson")[4:6]*100)))
printIQR = function(x) do.call(sprintf, c("%1.2g [IQR %1.2g \u2013 %1.2g]",as.list(quantile(x,c(0.5,0.25,0.75),na.rm = TRUE))))


summaryAnalysis = function(qualityAnalysisResult,lagAnalysisResult) {
  tmp4 = qualityAnalysisResult %>% ungroup() %>% select(model, comparison) %>% unnest(comparison)
  summ = tmp4 %>% group_by(model) %>% summarise(
    `average bias` = printIQR(residual),
    `average calibration` = printBinomCI(calibration),
    `average critical threshold` = printBinomCI(critical_threshold),
    `CRPS` = printIQR(crps),
    `quantile deviation` = sprintf("%1.3g", sqrt(sum((quantile.actual-0.5)^2,na.rm = TRUE)/(length(na.omit(quantile.actual))-1))-sqrt(1/12))
  )
  lag = lagAnalysisResult$modelLag %>% select(model,`estimate delay` = label)
  summ %>% inner_join(lag, by="model") %>% pivot_longer(cols = -model, values_to = "value",names_to="metric") %>% pivot_wider(values_from = "value", names_from="model")
}


plotRollingQuantiles = function(data, quantileVar, orderVar = "date", colours=c("red","blue"), window=14) {
  quantileVar = ensym(quantileVar)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  
  overallLag = data %>% summarise(median = quantile(!!quantileVar,0.5,na.rm = TRUE))
  
  summaryLag = data %>% arrange(!!orderVar) %>% 
    mutate(rollingQuant = slider::slide_index(!!quantileVar, !!orderVar, .before=window,.after=window,.f = ~ enframe(quantile(.x, c(0.025,0.25,0.5,0.75,0.975), na.rm=TRUE)))
    ) %>%
    select(!!!grps, !!orderVar,rollingQuant) %>%
    unnest(rollingQuant) %>%
    distinct()
  
  list(
    geom_line(data = summaryLag %>% filter(name %in% c("2.5%","97.5%")), mapping=aes(x=!!orderVar,y=value, group=name), linetype="dotted", colour=colours[1] ),
    geom_line(data = summaryLag %>% filter(name %in% c("25%","75%")), mapping=aes(x=!!orderVar,y=value, group=name), linetype="dashed", colour=colours[1] ),
    geom_line(data = summaryLag %>% filter(name %in% c("50%")), mapping=aes(x=!!orderVar,y=value, group=name), linetype="solid", colour=colours[1] ),
    geom_hline(data = overallLag, mapping = aes(yintercept = median), colour=colours[2])
  )
}

plotRollingDeciles = function(data, decileVar, orderVar = "date", colours="grey50", window=14) {
  decileVar = ensym(decileVar)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  
  summaryLag = data %>% arrange(!!orderVar) %>% 
    mutate(rollingQuant = slider::slide_index(!!decileVar, !!orderVar, .before=window,.after=window,.f = ~ enframe(quantile(.x, seq(0.1,0.9,length.out = 5), na.rm=TRUE)))
    ) %>%
    select(!!!grps, !!orderVar,rollingQuant) %>%
    unnest(rollingQuant) %>%
    distinct()
  
  list(
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=value, group=name), colour=colours, alpha=0.5)
  )
}

plotRollingProportion = function(data, binomialExpr, orderVar = "date", colours=c("red","blue"), window=14) {
  binomialExpr = enexpr(binomialExpr)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  
  overallLag = data %>% summarise(p = sum(!!binomialExpr,na.rm = TRUE)/length(na.omit(!!binomialExpr)))
  
  summaryLag = data %>% arrange(!!orderVar) %>% 
    mutate(
      x = slider::slide_index_dbl(!!binomialExpr, !!orderVar, .before=window,.after=window,.f = ~ sum(.x,na.rm = TRUE)),
      n = slider::slide_index_dbl(!!binomialExpr, !!orderVar, .before=window,.after=window,.f = ~ length(na.omit(.x))),
    ) %>%
    select(!!!grps, !!orderVar, x, n) %>%
    distinct() 
  
  # browser()
  
  summaryLag = summaryLag %>%
    filter(!is.na(x) & !is.na(n)) %>%
    mutate(binom::binom.confint(x,n,methods="wilson")) %>%
    rename(lower.0.025 = lower,upper.0.975 = upper) %>% select(-mean) %>%
    mutate(binom::binom.confint(x,n,methods="wilson",conf.level=0.5)) %>%
    rename(lower.0.25 = lower,upper.0.75 = upper)
  
  list(
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=lower.0.025), linetype="dotted", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=upper.0.975), linetype="dotted", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=lower.0.25), linetype="dashed", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=upper.0.75), linetype="dashed", colour=colours[1] ),
    geom_line(data = summaryLag, mapping=aes(x=!!orderVar,y=mean), linetype="solid", colour=colours[1] ),
    geom_hline(data = overallLag, mapping = aes(yintercept = p), colour=colours[2])
  )
}


plotRollingDensity = function(data, observationVar, n=50, orderVar = "date", colours=c("red","blue"), window=14) {
  observationVar = ensym(observationVar)
  orderVar = ensym(orderVar)
  grps = data %>% groups()
  obs = data %>% pull(!!observationVar)
  breaks = seq(min(obs,na.rm = TRUE),max(obs,na.rm = TRUE),length.out = n+1)
  miBreak = na.omit(breaks+lag(breaks))/2
  summaryLag = data %>% select(!!!grps,!!orderVar,!!observationVar) %>% nest(observations = !!observationVar) %>% arrange(!!orderVar) 
  summaryLag = summaryLag %>% mutate(bins = slider::slide(observations, .before=14, .after=14, .f = ~bind_rows(.x) %>% 
                                                            mutate(bin = cut(!!observationVar,breaks = breaks,labels = miBreak)) %>% 
                                                            group_by(bin) %>% 
                                                            summarise(count=n()) %>% 
                                                            mutate(prop = count/sum(count)) %>% 
                                                            ungroup() %>% 
                                                            tidyr::complete(bin, fill=list(count=0,prop=0)) ))
  summaryLag = summaryLag %>% select(!!!grps,!!orderVar,bins) %>% unnest(bins)
  
  list(
    geom_tile(data = summaryLag,  mapping=aes(x=!!orderVar, y=as.numeric(as.character(bin)), fill=prop))
  )
}