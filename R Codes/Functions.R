# Functions.
#=========================================================================
## Function to do the linear regression analysis ##
LMbyProtein<- function(data,dependent, independent){
  fit<-lm(y ~ x - 1,data)
  t<-summary(fit)
  
  data$predicted <- predict(fit)   # Save the predicted values
  data$residuals <- residuals(fit) # Save the residual values
  
  # Residuals graphic
  # Use the residuals to make an aesthetic adjustment 
  # (e.g. red colour when residual in very high) to highlight points 
  # which are poorly predicted by the model.
  Residuals<-abs(data$residuals)
  res_plot<-ggplot(data, aes(x = x, y = y)) +
    geom_smooth(method = "lm",formula = y~x-1, se = FALSE, color = "lightgrey") +
    geom_segment(aes(xend = x, yend = predicted), alpha = .2) +
    geom_point(aes(color = Residuals)) + # size also mapped
    geom_point(aes(y = predicted),shape = 1)+
    scale_color_continuous(low = "black", high = "red") +
    guides(color = guide_colorbar())+
    scale_x_continuous(breaks = seq(0,max(data$x)+0.1,by = 0.1))+
    scale_y_continuous(breaks = seq(0,max(data$y)+0.1,by = 0.1))+
    geom_label(x = min(data$x), hjust =0, y = max(data$y)-0.05,
               label = paste("RSE=", abbreviate(as.character(t$sigma),5),"\n",
                             "R-squared=", abbreviate(as.character(t$r.squared),5)))+
    ylab(dependent)+xlab(independent)+scale_fill_continuous(guide = guide_legend()) +
    theme(legend.key.width = unit(2.6, 'lines'), legend.position="bottom",
          axis.text.y =element_text(size=15),
          axis.text.x =element_text(size=15))
  
  # Coeff
  coef<-t$coefficients
  # Residuals
  residual<-t$residuals
  
  list(Res_plot=res_plot,coef=coef,res=residual)
}
#=========================================================================
# Remove outliers in data frame
remove_outliers <- function(data, na.rm = TRUE, ...) {
  x<-data$ratioOther
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE)
  H <- 1.5 * IQR(x, na.rm = na.rm)
  y <- x
  y[x < (qnt[1] - H)] <- NA
  y[x > (qnt[2] + H)] <- NA
  data$ratioOther<-y
  na.omit(data,cols="ratioOther")
}
#=========================================================================