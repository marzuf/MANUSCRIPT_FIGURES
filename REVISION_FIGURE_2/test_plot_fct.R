test_mat <- read.delim("test_matrix.txt", header=F)
test_mat[,2] <- 5
test_mat[3,] <- 10
test_mat[,5] <- 12

resol = 20000
test_mat2 <- test_mat
colnames(test_mat2) <- seq(from=20000, by=resol, length.out=ncol(test_mat2))
rownames(test_mat2) <- seq(from=20000, by=resol, length.out=nrow(test_mat2))

source("my_plot_matrix_v2.R")


my_plot_matrix(mat = test_mat2, 
               tad_coord = c(20001,220000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat2, 
               tad_coord = c(40001,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat2, 
               tad_coord = c(40001,120000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 2),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat2, 
               tad_coord = c(40001,100000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 3),
               main="", checkSim = F) 


test_mat <- read.delim("test_matrix.txt", header=F)
test_mat[,2] <- 5
test_mat[3,] <- 10
test_mat[,5] <- 12

resol = 20000

colnames(test_mat) <- seq(from=0, by=resol, length.out=ncol(test_mat))
rownames(test_mat) <- seq(from=0, by=resol, length.out=nrow(test_mat))

source("my_plot_matrix_v2.R")



my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,180000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 



my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,100000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 2),
               main="", checkSim = F) 







my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(60001,100000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(60001,100000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 3),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0,1),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 0),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 0),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,200000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 0),
               plotWithLeg=T,
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,40000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 1),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(0, 1),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,
               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 





source("my_plot_matrix.R")

my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,40000),
transformation = identity,               resolution = resol, 
               bins_around=c(0, 1),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 

my_plot_matrix_withLeg(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,               resolution = resol, 
               bins_around=c(0, 0),
               main="", checkSim = F) 

my_plot_matrix(mat = test_mat, 
               tad_coord = c(20001,60000),
               transformation = identity,               resolution = resol, 
               bins_around=c(1, 0),
               main="", checkSim = F) 



my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,60000),
               transformation = identity,               resolution = resol, 
               bins_around=c(1,0),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,60000),
               transformation = identity,               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,60000),
               transformation = identity,               resolution = resol, 
               bins_around=c(1, 1),
               main="", checkSim = F) 


my_plot_matrix(mat = test_mat, 
               tad_coord = c(1,200000),
               transformation = identity,               resolution = resol, 
               bins_around=c(0, 0),
               main="",  checkSim = F) 
