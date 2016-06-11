# An example of simulation studies for 2 popupations. 

# uniform distribution on p0
# estimate mf

source("funs_final.R")
library('ape')
library('ade4')
set.seed(1)

# number of simulation runs.
nsim = 10

# set p0 at the root of the tree to be uniform, with point mass at 0 and 1.
p0 = c(runif(80, 0, 1), rep(c(0, 1), each = 10))

# set branch length.
s2q1 = 0.1
s2q2 = 2*s2q1

# starting value for optimization.
sigmastart = c(0.1, 0.1)

# tolerance for optimization to stop.
steptol = 1.e-6

# record running time.
t = proc.time()

cat("s2q1 = ", s2q1, "; s2q2 = ", s2q2, "\n")
cat("nsim = ", nsim, "\n")

# run simulation and estimation.
res = sim2pop1(p0, nsim, s2q1, s2q2, sigmastart, zero40, weight40, steptol)
t = proc.time() - t
cat("\n", t, "\n")

# output estimated branch lengths. 
write(res, file = "./result/res2pop-0.1.txt", ncolumns = nsim)   

cat("mean = ", apply(res, 2, mean), "\n")
cat("var = ", apply(res, 2, var), "\n")

# second 2*100, L = 1000
p0 = c(runif(800, 0, 1), rep(c(0, 1), each = 100))

t = proc.time()
res = sim2pop1(p0, nsim, s2q1, s2q2, sigmastart, zero40, weight40, steptol)
t = proc.time() - t
cat("\n", t, "\n")

write(res, file = "./result/res2pop-0.1.txt", ncolumns = nsim,  append = TRUE) 

cat("mean = ", apply(res, 2, mean), "\n")
cat("var = ", apply(res, 2, var), "\n")

# third 2*100, L = 10000
p0 = c(runif(8000, 0, 1), rep(c(0, 1), each = 1000))

t = proc.time()
res = sim2pop1(p0, nsim, s2q1, s2q2, sigmastart, zero40, weight40, steptol)
t = proc.time() - t
cat("\n", t, "\n")

write(res, file = "./result/res2pop-0.1.txt", ncolumns = nsim,  append = TRUE) 

cat("mean = ", apply(res, 2, mean), "\n")
cat("var = ", apply(res, 2, var), "\n")

#############################################
# Simulation studies for n populations.
# African populations.

# true tree, from data analysis.
treeTrue = "(Mozabite:0.09783,(Mandenka:0.03682,(Yoruba:0.03171,((BiakaPygmy:0.03917,(MbutiPygmy:0.08769,San:0.46633):0.03809):0.02788,(BantuKenya:0.058,BantuSouthAfrica:0.09737):0.00677):0.02108):0.00226):0.09003);"
treeTrue = read.tree(text = treeTrue)

nsnp = 2000
steptol = 1.e-6

# record running time.
t = proc.time()
set.seed(1)


write.tree(treeSimRes$res, "./result/treeAfSim.txt")

write.tree(treeSimRes$resnj, "./result/treeAfSimUnroot.txt")

t = proc.time() - t

cat('\n', t, '\n')

# calculate consensus tree.
afCon = consensus(treeSimRes$res, p = 0.5)

# plot consensus tree.
plot(afCon)