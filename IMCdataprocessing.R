## To read file
library(EBImage)
library(dplyr)
library(ggplot2)

file.na <- readline(prompt = "Type the name of the .txt file, including extension: ")

main_df <- read.delim(file.na)

dimensions <- list( xnrow = max(main_df$X) + 1, yncol = max(main_df$Y) + 1)

print(paste("Image dimensions are: x", dimensions$xnrow,"um   y", dimensions$yncol, "um"  ))

marker_names <- data.frame(names(main_df)[7:length(names(main_df))])
index <- as.integer(1:(length(main_df)-6))

marker_names <- cbind(index, marker_names)
names(marker_names)[2] <- "marker"

nam <- as.character(marker_names$marker)
print(marker_names)

nuclei_marker <- readline(prompt = "Please type the number of the marker corresponding to the nuclei: ")

class(nuclei_marker) <- "integer"

print("Rescaling...")

rescaled_df <- sapply((7:ncol(main_df)), function(x) scales::rescale(main_df[,x], to = c(0,1)))

main_test <- cbind(main_df[,1:6], rescaled_df)
names(main_test) <- names(main_df)

for(i in 1:nrow(marker_names))
  base::assign(nam[i], value = matrix(main_test[,(6+i)], nrow = dimensions$xnrow, ncol = dimensions$yncol))

display_after_cutoff <- function(image_touse = nuclei_marker, LPD = .01, UPD = .99){
  x <- get(names(main_test)[image_touse+6])
  LVD <- quantile(x, probs = LPD)
  x[ x < LVD] <- LVD
  UVD <- quantile(x, probs = UPD)
  x[ x > UVD] <- UVD
  hist(x, breaks = 50)
  EBImage::display(x)
  return(x)
}


nuclear_segmentation <- function(x, thresh_w = 7, thresh_h = 7, xoffset = 0.0001, brushsize = 5, xshape = 'disc'){
  #nmask = watershed(distmap(x), 2)
  
  nmask = thresh(x, w=thresh_w,h=thresh_h, offset = xoffset)
  
  nmask = opening(nmask, makeBrush(brushsize, shape = xshape))
  
  nmask = fillHull(nmask)
  
  nmask = bwlabel(nmask)
  display(colorLabels(nmask))
  return(nmask)
}


cell_perimeter <- function(x, brush_size = a, xshape = 'diamond' ){
  kern = makeBrush(5, shape = xshape)
  dilated = dilate(x, kern)
  display(colorLabels(dilated))
  x <- as.vector(t(dilated))
  return(x)
}

reduce_data <- function(x, data = main_test ){
  test <- mutate(main_test, segments = x)
  test <- filter(test, segments != 0)
  median_values <- aggregate(test[,7:(length(names(main_df))-1)], list(test$segments), FUN = median)
  
  names(median_values)[1] <- "event"
  return(median_values)#
}

generate_fcs <- function(dta){
  filenam <- readline(prompt = "Type file name: ")
  library(Biobase)
  library(flowCore)
  library(flowViz)
  
dta <- as.matrix(dta)
  
  # you need to prepare some metadata
  meta <- data.frame(name=dimnames(dta)[[2]],
                     desc=paste('this is column',dimnames(dta)[[2]],'from your CSV')
  )
  meta$range <- apply(apply(dta,2,range),2,diff)
  meta$minRange <- apply(dta,2,min)
  meta$maxRange <- apply(dta,2,max)
  
  #head(meta)
  # all these are required for the following steps to work
  
  # a flowFrame is the internal representation of a FCS file
  ff <- new("flowFrame",
            exprs=dta,
            parameters=AnnotatedDataFrame(meta)
  )
  
  # a simple plot to check that it worked
  #xyplot(A~B,ff)
  
  # now you can save it back to the filesystem
  write.FCS(ff,paste(filenam,'FCS',sep='.'))
}
