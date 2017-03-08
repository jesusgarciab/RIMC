## To read file
library(EBImage)
library(dplyr)
library(ggplot2)

#to prompt user to type the file name
file.na <- readline(prompt = "Type the name of the .txt file, including extension: ")
#reading the file
main_df <- read.delim(file.na)
#getting the image dimensions from X and Y columns
dimensions <- list( xnrow = max(main_df$X) + 1, yncol = max(main_df$Y) + 1)

#Just displlaying the dimensions of the image
print(paste("Image dimensions are: x", dimensions$xnrow,"um   y", dimensions$yncol, "um"  ))

#getting the names of the channels
marker_names <- data.frame(names(main_df)[7:length(names(main_df))])
#using an index so the person can easily "select" what channel to use next
index <- as.integer(1:(length(main_df)-6))
marker_names <- cbind(index, marker_names)
names(marker_names)[2] <- "marker"
nam <- as.character(marker_names$marker)
print(marker_names)

#asking for the intercalator/nuclei channel for use later on
nuclei_marker <- readline(prompt = "Please type the number of the marker corresponding to the nuclei: ")
class(nuclei_marker) <- "integer"

#not sure.. but I thought it would look cool
print("Rescaling...")

#rescaling to 0-1
rescaled_df <- sapply((7:ncol(main_df)), function(x) scales::rescale(main_df[,x], to = c(0,1)))
#getting the whole initial data.frame info but now with rescaled data
main_test <- cbind(main_df[,1:6], rescaled_df)
names(main_test) <- names(main_df)

#converting each column into matrices with the specified dimensions
for(i in 1:nrow(marker_names))
  base::assign(nam[i], value = matrix(main_test[,(6+i)], nrow = dimensions$xnrow, ncol = dimensions$yncol))

#Function to display the image of interest (default nuclei_marker) and cropping by upper and lower percentiles
display_after_cutoff <- function(image_touse = nuclei_marker, LPD = .01, UPD = .99){
  #get the variable name and assign to x
  x <- get(names(main_test)[image_touse+6])
  #to calculate lower value of detection
  LVD <- quantile(x, probs = LPD)
  x[ x < LVD] <- LVD
  #to calculate upper limit of detection
  UVD <- quantile(x, probs = UPD)
  x[ x > UVD] <- UVD
  #histogram of the image
  hist(x, breaks = 50)
  #display the image
  EBImage::display(x)
  return(x)
}

#function for nuclear segmentation (the default values are the ones that I came up with after a little tweaking)
nuclear_segmentation <- function(x, thresh_w = 7, thresh_h = 7, xoffset = 0.0001, brushsize = 5, xshape = 'disc'){
  #nmask = watershed(distmap(x), 2)
  #Copied from an EBimage example to come up with an image of only white/black pixels
  nmask = thresh(x, w=thresh_w,h=thresh_h, offset = xoffset)
  
  nmask = opening(nmask, makeBrush(brushsize, shape = xshape))
  
  nmask = fillHull(nmask)
  
  nmask = bwlabel(nmask)
  display(colorLabels(nmask))
  return(nmask)
}

#function to "dilate" nucleus. I don't know why i used the 5 pixel brush...
cell_perimeter <- function(x, brush_size = a, xshape = 'diamond' ){
  kern = makeBrush(5, shape = xshape)
  dilated = dilate(x, kern)
  display(colorLabels(dilated))
  x <- as.vector(t(dilated))
  return(x)
}
#function to get the median of the pixels within the area of each "cell area"
reduce_data <- function(x, data = main_test ){
  #add segments(cells) column
  test <- mutate(main_test, segments = x)
  #eliminating data of no cells
  test <- filter(test, segments != 0)
  #get meadian
  median_values <- aggregate(test[,7:(length(names(main_df))-1)], list(test$segments), FUN = median)
  
  names(median_values)[1] <- "event"
  return(median_values)#
}

#function to generate fcs file
generate_fcs <- function(dta){
  filenam <- readline(prompt = "Type file name: ")
  library(Biobase)
  library(flowCore)
  library(flowViz)
  
dta <- as.matrix(dta)
  
  # preparing some metadata
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
  

  
  # now you can save it back to the filesystem
  write.FCS(ff,paste(filenam,'FCS',sep='.'))
}
