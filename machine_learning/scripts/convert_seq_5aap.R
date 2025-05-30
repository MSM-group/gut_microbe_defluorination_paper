convert_seq_5aap <- function (x) 
{
  tmp <- read_csv("data/5_aa_properties.csv") %>%
    dplyr::filter(complete.cases(.)) 
  tmpcut <- tmp %>% 
    select(-1) 
  tmpfin <- data.frame(sapply(tmpcut, as.numeric))
  
  aalist <- list()
  for(i in 1:nrow(tmpfin)) {
    aalist[[i]] <- as.numeric(dplyr::slice(tmpfin, i))
  }
  names(aalist) <- as.vector(tmp$AA_ABREV)
  
  gap <- rep(0, 5)
  
  z <- length(unlist(strsplit(x, split = NULL)))
  ans <- vector()
  for (i in 1:z) {
    if ((x[i] == "D") || (x[i] == "d")) {
      ans <- c(ans, aalist[["D"]])
    }
    if ((x[i] == "T") || (x[i] == "t")) {
      ans <- c(ans, aalist[["T"]])
    }
    if ((x[i] == "E") || (x[i] == "e")) {
      ans <- c(ans, aalist[["E"]])
    }
    if ((x[i] == "C") || (x[i] == "c")) {
      ans <- c(ans, aalist[["C"]])
    }
    if ((x[i] == "M") || (x[i] == "m")) {
      ans <- c(ans, aalist[["M"]])
    }
    if ((x[i] == "Y") || (x[i] == "y")) {
      ans <- c(ans, aalist[["Y"]])
    }
    if ((x[i] == "K") || (x[i] == "k")) {
      ans <- c(ans, aalist[["K"]])
    }
    if ((x[i] == "R") || (x[i] == "r")) {
      ans <- c(ans, aalist[["R"]])
    }
    if ((x[i] == "S") || (x[i] == "s")) {
      ans <- c(ans, aalist[["S"]])
    }
    if ((x[i] == "Q") || (x[i] == "q")) {
      ans <- c(ans, aalist[["Q"]])
    }
    if ((x[i] == "F") || (x[i] == "f")) {
      ans <- c(ans, aalist[["F"]])
    }
    if ((x[i] == "P") || (x[i] == "p")) {
      ans <- c(ans, aalist[["P"]])
    }
    if ((x[i] == "W") || (x[i] == "w")) {
      ans <- c(ans, aalist[["W"]])
    }
    if ((x[i] == "N") || (x[i] == "n")) {
      ans <- c(ans, aalist[["N"]])
    }
    if ((x[i] == "G") || (x[i] == "g")) {
      ans <- c(ans, aalist[["G"]])
    }
    if ((x[i] == "V") || (x[i] == "v")) {
      ans <- c(ans, aalist[["V"]])
    }
    if ((x[i] == "I") || (x[i] == "i")) {
      ans <- c(ans, aalist[["I"]])
    }
    if ((x[i] == "L") || (x[i] == "l")) {
      ans <- c(ans, aalist[["L"]])
    }
    if ((x[i] == "A") || (x[i] == "a")) {
      ans <- c(ans, aalist[["A"]])
    }
    if ((x[i] == "H") || (x[i] == "h")) {
      ans <- c(ans, aalist[["H"]])
    }
    if (x[i] == "-") {
      ans <- c(ans, gap)
    }
    if (x[i] == ".") {
      ans <- c(ans, gap)
    }
    if (x[i] == "X") {
      ans <- c(ans, gap)
    }
    if (!length(ans) == 5 * i) {
      print(c("Error for position ", i, "in ", x))
      return("ERROR")
    }
  }
  return(ans)
}
