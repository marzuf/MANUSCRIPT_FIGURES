my_func <- function(x=5) {

z <- foreach(i=1:x, .combine='c') %dopar% {
	return(i)
}
return(z)
}
