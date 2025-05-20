int_score <- function(alpha, truth, L, U){
  return(
    (U - L) + 2/alpha*(truth < L)*(L - truth) + 2/alpha*(truth > U)*(truth - U)
  )
}

