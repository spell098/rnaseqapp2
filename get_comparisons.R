#' Take all the name and export all possible comparisons
#' @author Simon J Pelletier
#' @param names Names of all groups compared
#' @return A vector of all the possible comparisons in the format names1//name2 (VERIFY!!)
#' @keywords comparisons
#' @export
get_comparisons = function(names){
  comparisons=NULL
  x <- strsplit(as.character(names[i]),"_")[[1]]
  for(i in 1:length(x)){
    if(!is.na(as.numeric(x[i]))){
      l <- i-1
    } 
  }
  conditions = matrix(nrow=length(names),ncol = l)
  for (i in 1:length(names)){
    conditions[i,] = strsplit(names[i],"_")[[1]][1:l]
  }
  for (j in 1:ncol(conditions)){
    tmp=unique(conditions[,j])
    comparisons[j] = paste(tmp,collapse="\\")
  }
  return(comparisons)
}
