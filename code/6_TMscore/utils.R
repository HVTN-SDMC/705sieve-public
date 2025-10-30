library(ggrepel)
library(EnvStats)

volcano_plot <- function(df, show_labels = TRUE){
  p_breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  
  p <- ggplot(df, aes(x = -markCoef, y = -log10(p.unadj))) +
    geom_point(aes(color = p.FDR), size = 1.5, alpha = 1,
               shape = 21, stroke = 2, fill = "white") +
    scale_x_continuous(breaks = -5:5) +
    scale_y_continuous(breaks = -log10(p_breaks), labels = p_breaks, 
                       minor_breaks = NULL, limits = c(0, -log10(0.0017))) +
    labs(x = "Difference in Log Hazard Ratio (Vaccine/Placebo)\nfor Min vs. Max TM Score",
         y = "Unadjusted Double One-Sided Sieve P-value",
         title = paste0("Protomer: ", unique(df$protomer), "; lower bound: ", 
                        substring(unique(df$lb), first = 3))) +
    # scale_color_gradient(name = "Q-value",breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    #                      labels = c("0", "0.2", "0.4", "0.6","0.8", "1"),
    #                      limits = c(0,1), high = "green", low = "blue")+
    scale_color_viridis_c(name = "Q-value",breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                         labels = c("0", "0.2", "0.4", "0.6","0.8", "1"),
                         limits = c(0,1), option = "H")+
    annotate("segment", x = -0.2, xend = -3, y = -log10(0.003), yend = -log10(0.003), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("segment", x = 0.2, xend = 3, y = -log10(0.003), yend = -log10(0.003), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("text", x = -1.7, y = -log10(0.002), label = "Greater protection against\nlower TM score", size = 4)+
    annotate("text", x = 1.6, y = -log10(0.002), label = "Greater protection against\nhigher TM score", size = 4)+
    
    theme_bw()+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(hjust = 0, vjust = 2, size = 14),
          axis.title.x = element_text(size = 14, vjust = -2),
          axis.title.y = element_text(size = 14, vjust = 2),
          axis.text.x = element_text(size = 14, colour = "black"),
          axis.text.y = element_text(size = 14, colour = "black"),
          legend.text=element_text(size=14),
          legend.title = element_text(size=14, vjust =3 ),
          legend.key.size = unit(1.5, "line"))
  
  if (show_labels){
    p <- p + 
      geom_text_repel(aes(label = lab), size = 3, force = 15, seed = 0,
                      #data = filter(df, pHRconstancy < 0.007), 
                      data = filter(df, p.FDR <= 0.2),
                      max.overlaps = Inf)
  }
  
  return(p)
}

format.p <- function(p, ndigits=2){
  pp <- NULL
  for(i in 1:length(p)){
    if(is.na(p[i])){
      pp[i] <- "--"
    }else{
      if(p[i]<0.001){pp[i] <- " < 0.001"}
      else if (p[i]==1){pp[i] <- "= 1"}
      else{pp[i] <-paste0(" = ",as.character(format(as.numeric(p[i]),digits=ndigits,nsmall=ndigits)))
      if(pp[i]== " = 1.00") {pp[i]= " = 1"}}
    }
  }
  return (pp)
}

volcano_plot_beta <- function(df, show_labels = TRUE){
  p_breaks <- c(0.001, 0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1)
  
  p <- ggplot(df, aes(x = beta, y = -log10(p.unadj))) +
    geom_point(aes(color = p.FDR), size = 3, alpha = 1,
               shape = 21, stroke = 2, fill = "white") +
    scale_x_continuous(breaks = c(-0.1, -0.05, 0, 0.05, 0.1), limits = c(-0.11, 0.11),
                       label = c("-0.1", "-0.05", "0", "0.05", "0.1")) +
    scale_y_continuous(breaks = -log10(p_breaks), labels = p_breaks, 
                       minor_breaks = NULL, limits = c(0, -log10(0.008))) +
    labs(x = "Mean Difference in Standarized TM Score\n(Vaccine vs. Placebo)",
         y = "Unadjusted Two-Sided P-value",
         title = paste0("Protomer: A", "; lower bound: 6")) +
    # scale_color_gradient(name = "Q-value",breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
    #                      labels = c("0", "0.2", "0.4", "0.6","0.8", "1"),
    #                      limits = c(0,1), high = "green", low = "blue")+
    scale_color_viridis_c(name = "Q-value",breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1),
                          labels = c("0", "0.2", "0.4", "0.6","0.8", "1"),
                          limits = c(0,1), option = "H")+
    annotate("segment", x = -0.02, xend = -0.1, y = -log10(0.017), yend = -log10(0.017), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("segment", x = 0.02, xend = 0.1, y = -log10(0.017), yend = -log10(0.017), linewidth = 1, arrow = arrow(length = unit(0.3, "cm")))+
    annotate("text", x = -0.06, y = -log10(0.011), label = "Placebo Structures Closer to\nVaccine Insert Structure", size = 4.7)+
    annotate("text", x = 0.06, y = -log10(0.011), label = "Vaccine Structures Closer to\nVaccine Insert Structure", size = 4.7)+

    theme_bw()+
    theme(plot.margin = unit(c(1,1,1,1), "cm"),
          plot.title = element_text(hjust = 0, vjust = 2, size = 17),
          axis.title.x = element_text(size = 17, vjust = -2),
          axis.title.y = element_text(size = 17, vjust = 2),
          axis.text.x = element_text(size = 17, colour = "black"),
          axis.text.y = element_text(size = 17, colour = "black"),
          legend.text=element_text(size=17),
          legend.title = element_text(size=17, vjust =3 ),
          legend.key.size = unit(1.5, "line"))
  
  if (show_labels){
    p <- p + 
      geom_text_repel(aes(label = lab), size = 3, force = 15, seed = 0,
                      #data = filter(df, pHRconstancy < 0.007), 
                      data = filter(df, p.FDR <= 0.2),
                      max.overlaps = Inf)
  }
  
  return(p)
}
