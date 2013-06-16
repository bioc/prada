## Functions for dealing with alphanumeric identifiers for larger well plates
## There may be one or two letters in the string
getAlphaNumeric = function(horizontal, vertical) {
    if (any(horizontal>702)) stop(sprintf("Indices of 'horizontal' well must not exceed %d", 26*27))
    if (any(vertical>99)) stop(sprintf("Indices of 'horizontal' well must not exceed %d.", 99))
    
    alpha1 <- c("", LETTERS) [(horizontal - 1)%/%length(LETTERS) + 1]
    alpha2 <- LETTERS[(horizontal-1) %% length(LETTERS) +1]
    id.num <- sprintf('%02d', vertical)
    id.alpha <- paste(alpha1, alpha2, sep='')
    id <- paste(id.alpha, id.num, sep='')
    return(list(id=id, id.alpha=id.alpha, id.num=id.num))
}
