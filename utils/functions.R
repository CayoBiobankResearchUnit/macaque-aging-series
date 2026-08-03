jaccard<- function(set1, set2) {
  intersection <- length(intersect(set1, set2))
  union <- length(unique(c(set1, set2)))
  if (union == 0) return(0)  # Handle edge case
  return(intersection / union)
}