#options(warn=2)
#options(warn=1)

## Lag analysis ----

#' Title
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
  triangular = tibble(name = triangularSim$name, ts = list(triangularSim$ts), infectivityProfile = list(triangularSim$infectivityProfile$yMatrix))
  
  periodicEstimates = triangular %>% inner_join(estimators, by=character()) %>% mutate(
    estimate = pmap(.l=list(fn = estimFn, ts = ts, infectivityProfile = infectivityProfile), .f = function(fn,ts,infectivityProfile) {
      #browser()
      fn(ts,infectivityProfile)
    })
  )

  periodicRt = periodicEstimates %>% select(model,estimate) %>% unnest(estimate)
  periodicTs = periodicEstimates %>% select(model,ts) %>% unnest(ts)

  offset = tibble(shift = c(0:28)) %>% inner_join(
    periodicTs %>% select(theoretical.date = date, Rt.actual), by=character()) %>% 
    mutate(
      estimate.date = theoretical.date+shift
    )
  
  offsetPeriodicRt = periodicRt %>% 
    inner_join(offset, by=c("date"="estimate.date"))
  
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
    estimates = periodicRt %>% inner_join(periodicTs, by=c("model","date","subgroup")),
    rmseByOffset = offsetSummary,
    modelLag = offsetSummarySummary
  )
  
  class(out) = "lag_analysis_result"
  
  return(out)
}

## Lag analysis visualisation ----

lagTimeseriesPlot = function(lagAnalysisResult, estimation = "Rt", ncol=1) {
  estimates = lagAnalysisResult$estimates
  trueVar = as.symbol(paste0(estimation, ".actual"))
  estimateVar = as.symbol(paste0(estimation, ".Quantile.0.5"))
  lowVar = as.symbol(paste0(estimation, ".Quantile.0.025"))
  highVar = as.symbol(paste0(estimation, ".Quantile.0.975"))
  
  defaultAes = aes(x=date, y=Rt.Mean)
  dots = rlang::list2(...)
  combinedAes = modifyList(defaultAes,mapping)
  
  suppressWarnings(
    ggplot(estimates %>% filter(!is.nan(!!trueVar)), mapping = aes(x=date))+
      geom_ribbon(aes(ymin = !!lowVar, ymax = !!highVar), alpha=0.2,colour = NA,fill="grey20")+
      geom_line(aes(y=!!estimateVar),colour="black")+
      guides(colour="none")+
      geom_line(mapping = aes(x=date,y=!!trueVar), inherit.aes = FALSE, colour="red")+
      xlab(NULL)+
      ylab(paste0(estimation," estimate"))+
      facet_wrap(vars(model),ncol=ncol)
  )
}

lagOffsetPlot = function(rmseByOffset, modelLag, ncol=1) {
  ggplot(rmseByOffset, aes(x=shift,y=rmse))+
    geom_point()+
    geom_line()+
    geom_vline(data = modelLag, mapping=aes(xintercept = medianLag), colour="blue")+
    geom_label(data = modelLag, mapping=aes(label=label,x = medianLag), y=Inf, fill=NA, label.size=0, label.padding = unit(1, "lines"),hjust=0, vjust="inward", colour="blue",size=standardPrintOutput::labelInPoints(8))+
    ylab(latex2exp::TeX("$R_t$ RMSE"))+xlab("shift (days)")+
    facet_wrap(vars(model),ncol=ncol)
}

lagPlot = function(lagAnalysis) {
  p1 = lagTimeseriesPlot(lagAnalysis$estimates, ncol=1)
  p2 = lagOffsetPlot(lagAnalysis$rmseByOffset, lagAnalysis$modelLag, ncol=1)
  p3 = p1+p2+plot_layout(nrow=1,widths = c(2,1))+plot_annotation(tag_levels = "A")
  return(p3)
}


# reconstruct a CDF on a regular grid from the quantiles
inferQuantiles = function(rtEstimate, prefix="Rt.") {
  nm = function(s) as.symbol(paste0(prefix,s))
  quantilePrefix = paste0(prefix,"Quantile.")
  rtEstimate %>% 
    mutate(
      median = !!nm("Quantile.0.5"),
      lo = !!nm("Quantile.0.025"),
      hi = !!nm("Quantile.0.975"),
      iqr = !!nm("Quantile.0.75")-!!nm("Quantile.0.25")
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

calculateMetrics = function(quantileEstimate, originalTs, medianLag, trueVar = "Rt.actual", criticalValue = 1) {
  trueVar = ensym(trueVar)
  originalTs %>% 
    select(subgroup,date, actual = !!trueVar) %>%
    left_join( # here is the point where we truncate the data set as a result of shifting the estiamtes by the lag. If this were a left join
      quantileEstimate %>% 
        mutate(
          estimate.date = date,
          date = estimate.date-round(medianLag)),
      by=c("subgroup","date")
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
    ) %>% mutate(
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

qualityAnalysis = function(estimators, modelLag, estimatePrefix = "Rt.", trueVar = "Rt.actual", criticalValue = 1, synthetic = jepidemic::validationDataset) {
  trueVar = ensym(trueVar)
  
  syntheticEstimates = synthetic %>% inner_join(estimators, by=character()) %>% mutate(
    estimate = pmap(.l=list(fn = estimFn, ts = ts, infectivityProfile = infectivityProfile), .f = function(fn,ts,infectivityProfile) {
      fn(ts,infectivityProfile)
    })
  )
  
  syntheticEstimates2 = syntheticEstimates %>% 
    mutate(
      quantileEstimate = map(estimate, inferQuantiles, prefix=estimatePrefix)
    ) %>%
    inner_join(modelLag, by="model") %>%
    mutate(
      comparison = pmap(.l=list(quantileEstimate, ts, medianLag), .f=calculateMetrics, trueVar=!!trueVar, criticalValue = criticalValue)
    ) %>%
    select(-quantileEstimate)
  
  
}

## Validation plots ----

estimationExamplePlot = function(qualityAnalysisResult, simFilterExpr = {weekendEffect == 0.1 & seed == 100}, modelFilterExpr = TRUE, subgroups = 1) {
  simFilterExpr = enexpr(simFilterExpr)
  modelFilterExpr = enexpr(modelFilterExpr)
  
  tmp = qualityAnalysisResult %>% 
    filter(!!simFilterExpr) %>%
    filter(!!modelFilterExpr) %>%
    mutate(label = paste0(smoothLabel,"; ",seedLabel,"; ",weekendLabel))
  
  tmp2 = tmp %>% ungroup() %>% select(label,ts) %>% distinct() %>% unnest(ts) %>% filter(subgroup %in% subgroups) %>% mutate(row="simulation")
  tmp3 = tmp %>% ungroup() %>% select(model,label,estimate) %>% unnest(estimate) %>% filter(subgroup %in% subgroups)
  
  p1 = incidencePlot(tmp2)+geom_line(aes(y=lambda_t),colour="red")+facet_grid(rows = vars(row), cols=vars(label))
  
  p2 = rtRibbon(tmp3,ylim = c(0.6,1.4), aes(group=subgroup))+
    geom_line(data=tmp2,mapping=aes(x=date,y=Rt.actual),colour="red")+
    ylab(latex2exp::TeX("$R_t$ estimate"))+xlab(NULL) +facet_grid(rows = vars(model),cols = vars(label))
  
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
  p1 = ggplot(tmp3 %>% filter(subgroup==1), aes(x=date))+
    geom_point(aes(x=date,y=value),size=0.1)+ 
    geom_line(aes(y=lambda_t), colour="red")+
    facet_wrap(vars(label),nrow = 1)+ylab(latex2exp::TeX("model cases"))+scale_y_continuous(trans="log1p",breaks = ukcovidtools::breaks_log1p())
  
  p2 = ggplot(tmp4 %>% filter(subgroup==1), aes(x=date))+
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