

#Function written by Ruben Props
ram_contrast <- function(hyprs, comp1, comp2, plot = TRUE){
  # Extract spectra
  x <- hyprs@data$spc
  
  # Remove clusters column
  Samples <- factor(unlist(hyprs@label))
  Samples <- Samples[Samples!="clusters"]; Samples <- droplevels(Samples)
  lev1 <- comp1
  lev2 <- comp2
  comp1 <- Samples %in% comp1
  comp2 <- Samples %in% comp2
  
  # Print table of samples
  cat(paste0("-----------------------------------------------------------------------------------------------------", 
             "\n \n"))
  cat(paste0("\t Your cells are distributed over these samples:\n\n "))
  print(table(Samples))
  cat(paste0("-----------------------------------------------------------------------------------------------------", 
             "\n \n"))
  max.total <- max(x)
  if (length(comp1) == 1 & length(comp2) == 1) 
    tmp <- (x[comp1, ] - x[comp2, ])
  if (length(comp1) == 1 & length(comp2) != 1) 
    tmp <- x[comp1, ] - colMeans(x[comp2, ])
  if (length(comp1) != 1 & length(comp2) == 1) 
    tmp <- colMeans(x[comp1, ]) - x[comp2, ]
  if (length(comp1) != 1 & length(comp2) != 1) 
    tmp <- colMeans(x[comp1, ]) - colMeans(x[comp2,])
  tmp <- tmp/max.total
  results <- data.frame(Density = tmp, Wavenumber = hyprs@wavelength)
  cat(paste0("\t Returning contrasts between mean spectra for ", sum(comp1),
             " cells of\n ", list(lev1)))
  cat(paste0("\n\t and ", sum(comp2) ," cells of\n ",
             list(lev2), "\n"))
  cat(paste0("-----------------------------------------------------------------------------------------------------", 
             "\n \n"))
  if(plot){
    # Plot
    v <- ggplot2::ggplot(results, ggplot2::aes(x=Wavenumber, y=Density))+
      geom_line()+
      labs(x=bquote("Wavenumber [cm"^-1*"]"),y="Density[A.U.]")+
      ggplot2::theme(axis.text=element_text(size=15),
                     axis.title=element_text(size=18))+
      #ggplot2::scale_fill_distiller(palette="RdBu", na.value="white") +
      ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                     panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      ggplot2::ylim(-0.38, 0.38)
    print(v)
  }
  return(results)
}


###############################################
#fancy y axis
fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}
