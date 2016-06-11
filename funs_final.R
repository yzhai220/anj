# Functions used in the ANJ project.
library("ape")
library("hash")
library("statmod")

# simNormAppr
# approximated truncated normal
# i.e., putting two mass atom at 0 and 1, and using sigma

simNormAppr = function(p, s2q)
{
  if (s2q < 0)
  {
    cat("Negative branch length sets to 0. \n")
    s2q = 0
  }
  sigma = sqrt(s2q)
  sd0 = sqrt(p * (1 - p))
  res = rnorm(length(p), mean = p, sd = sd0 * sigma)
  res[res > 1] = 1
  res[res < 0] = 0
  res
}


# denp0vec1
# density of p0, to be integrated, for any 0 < p0 < 1
# for vecotrs p1 only
# uniform prior over [0, 1]

denp0vec1 = function(p0vec, s2q1, p1vec)
{
  np0 = length(p0vec)
  np1 = length(p1vec)
  sigma1 = sqrt(s2q1)

  sdp0 = sqrt(p0vec * (1 - p0vec))
  sdp1 = rep(sdp0 * sigma1, each = np1)
  mup0 = rep(p0vec, each = np1)

  p1index0 = rep(p1vec == 0, np0)
  p1index1 = rep(p1vec == 1, np0)
  p1index2 = !(p1index0 | p1index1)

  p1den = rep(NA, np1 * np0)

  p1den[p1index0] = pnorm(0, mean = mup0[p1index0], sd = sdp1[p1index0])
  p1den[p1index1] = pnorm(1, mean = mup0[p1index1], sd = sdp1[p1index1],
                          lower.tail = F)
  p1den[p1index2] = dnorm(rep(p1vec, np0)[p1index2], mean = mup0[p1index2],
                          sd = sdp1[p1index2])

  p1den
}


# denp0vec
# density of p0, to be integrated, for any 0 < p0 < 1
# for vecotrs p1, p2
# uniform prior over [0, 1]

denp0vec = function(p0vec, s2q1, s2q2, p1vec, p2vec)
{
  np0 = length(p0vec)
  np1 = length(p1vec)

  p1den = denp0vec1(p0vec, s2q1, p1vec)
  p2den = denp0vec1(p0vec, s2q2, p2vec)

  den = p1den * p2den
  den = matrix(den, np0, np1, byrow = TRUE)
  den
}


# nllksigmavec
# negative log-likelihood using denp0vec

nllksigmavec = function(param, p1, p2, zero, weight, mf)
{
  np1 = length(p1)
  s2q1 = param[1]
  s2q2 = param[2]
  if (s2q1 <= 0 | s2q2 <= 0) { return(1.e10) }

  mn = 1 - mf

  denmatrix = denp0vec(zero/2 + 0.5, s2q1, s2q2, p1, p2) * (weight/2)
  tem = colSums(denmatrix)

  index0 = (p1 == 0 & p2 == 0)
  index1 = (p1 == 1 & p2 == 1)
  indexf = index0 | index1

  tem = tem * mn   # scale the continuous density
  tem[indexf] = tem[indexf] + mf/2   # scale the fixation probability

  nllk = -sum(log(tem))
  nllk
}


# nllksurf
# nllk surface for different sigma values

nllksurf = function(grid1, grid2, p1, p2, zero, weight)
{
  npt1 = length(grid1)
  npt2 = length(grid2)
  z = matrix(NA, npt1, npt2)

  for (i in 1:npt1)
  {
    for (j in 1:npt2)
      z[i, j] = nllksigmavec(c(grid1[i], grid2[j]), p1, p2, zero, weight)
  }
  z
}


# Gaussian quadrature points.

quad20 = gauss.quad(20, kind = "legendre")
zero20 = quad20$nodes
weight20 = quad20$weights

quad40 = gauss.quad(40, kind = "legendre")
zero40 = quad40$nodes
weight40 = quad40$weights

quad80 = gauss.quad(80, kind = "legendre")
zero80 = quad80$nodes
weight80 = quad80$weights


# symm
# index indicates those should x be changed to 1 - x
# alternatively change x to 1-x

symm = function(data)
{
  ndata = length(data)
  index = rep(c(T, F), ndata/2)

  data[index] = 1 - data[index]
  data
}


# distm
# function to calculate the distances of all pairs
# the distance of a pair is an array (s2q1, s2q2)
# for pi and pj, the distance s2q1 is stored at (i, j)
# and s2q2 is stored at (j, i)

# input: a n*l data matrix, with rows for pops, and cols for alleles
# output: a n*n distance matrix
# ith column of the matrix represents the distances from pi,
# to the ancestor of pi and another pop

distm = function(data, sigmastart = c(0.2, 0.2), zero = zero40, weight = weight40,
                 mf = 0.2)
{
  n = dim(data)[1]
  L = dim(data)[2]

  res = diag(NA, n)

  for (i in 1:(n-1))
  {
    cat(i, '\n')
    temp1 = data[i, ]
    for (j in (i+1):n)
    {
      temp2 = data[j, ]
      tem = nlm(nllksigmavec, sigmastart, p1 = temp1, p2 = temp2, zero = zero,
                weight = weight, mf = mf, print.level = 0,
                iterlim = 200, steptol = 1.e-3)
      if (tem$code != 1 && tem$code != 2)
      { cat('Warning: not converge for', i, 'and', j, 'Code:', tem$code, '\n')}
      res[j, i] = tem$est[1]
      res[i, j] = tem$est[2]
    }
  }
  res
}


# distm1
# use mf calculated from data, instead of fixed number
# naive estimate of mf

distm1 = function(data, sigmastart = c(0.2, 0.2), zero = zero40,
                  weight = weight40, steptol = 1.e-5)
{
  n = dim(data)[1]
  L = dim(data)[2]

  res = diag(NA, n)

  # Warm start, or different start values for different pairs.
  if (length(sigmastart) == n*n)
  {
    customerized_sigmastart = TRUE
  }
  else
  {
    customerized_sigmastart = FALSE
  }

  for (i in 1:(n-1))
  {
    cat(i, '\n')
    temp1 = data[i, ]
    temp1 = symm(temp1)
    for (j in (i+1):n)
    {
      temp2 = data[j, ]
      temp2 = symm(temp2)

      datat = (temp1 + temp2)/2
      index.fix0 = datat == 0
      index.fix1 = datat == 1
      m0 = sum(index.fix0) / L
      m1 = sum(index.fix1) / L
      mf = m0 + m1

      if (customerized_sigmastart == TRUE)
      {
        sigmastartij = c(sigmastart[j, i], sigmastart[i, j])
      }
      else
      {
        sigmastartij = sigmastart
      }
      tem = nlm(nllksigmavec, sigmastartij, p1 = temp1, p2 = temp2, zero = zero,
                weight = weight, mf = mf, print.level = 0,
                iterlim = 200, steptol = steptol)
#      if (tem$code != 1 && tem$code != 2)
#      { cat('Warning: not converge for', i, 'and', j, 'Code:', tem$code, '\n')}
      res[j, i] = tem$est[1]
      res[i, j] = tem$est[2]
    }
  }
  res
}


# tree0
# ANJ tree

tree0 = function(distmat)
{
  distmatn = distmat
  n = dim(distmatn)[1]
  name = 1:n
  nnodes = c(rep(1, n), rep(NA, n-1))

  topo = matrix(NA, n-1, 4)
  leng = matrix(NA, n-1, 2)

  for (i in 1:(n-1))
  {
    #    cat(i)    # monitor the appearance of bad trees, i.e., loops in the tree
    index = which.min(distmatn)
    temp = index / (n - i + 1)
    nj = ceiling(temp)
    ni = index - (nj-1) * (n - i + 1)
    indexn = which.min(distmatn[, ni])
    iter = 1
    while ((indexn != nj) && (iter <= (n - i + 1)))
    {
      nj = ni
      ni = indexn
      indexn = which.min(distmatn[, ni])
      iter = iter + 1         # control for bad trees
    }
    tempindex = c(ni, nj)
    namenew = n+i
    topo[i, 1:2] = name[tempindex]   # 1st & 2nd columns: nodes combined
    topo[i, 3] = namenew    # 3rd column: the number of the new node

    nnodes1 = nnodes[name[ni]]
    nnodes2 = nnodes[name[nj]]
    nnodes[n+i] = nnodes1 + nnodes2
    topo[i, 4] = nnodes[n+i]    # 4th column: the number of the leaves in this new node

    leng[i, ] = c(distmatn[nj, ni], distmatn[ni, nj])
    #    temleng = mean(distmatn[, ni] - distmatn[, nj], na.rm = TRUE)
    #    distmatn[, ni] - temleng
    #    resleng[i, 1] = c(distmatn[nj, ni], distmatn[ni, nj])

    meandist1 = (distmatn[, ni] + distmatn[, nj]
                 - distmatn[nj, ni] - distmatn[ni, nj]) / 2
    distmatn[, nj] = meandist1
    meandist2 = (distmatn[ni, ] + distmatn[nj, ]) / 2
    distmatn[nj, ] = meandist2
    name[nj] = namenew
    distmatn = distmatn[-ni, -ni]
    name = name[-ni]
#    distmatn = pmax(distmatn, 0)       # make all branch lengthes positive
  }
  list(topo = topo, leng = leng)
}


# tree1
# Neighbor-Joining (NJ) updating rule, in the asymmetric fashion.

tree1 = function(distmat)
{
  distmatn = distmat
  n = dim(distmatn)[1]
  name = 1:n
  nnodes = c(rep(1, n), rep(NA, n-1))

  topo = matrix(NA, n-1, 4)
  leng = matrix(NA, n-1, 2)

  for (i in 1:(n-1))
  {
    #    cat(i)    # monitor the appearance of bad trees, i.e., loops in the tree
    index = which.min(distmatn)
    temp = index / (n - i + 1)
    nj = ceiling(temp)
    ni = index - (nj-1) * (n - i + 1)
    indexn = which.min(distmatn[, ni])
    iter = 1
    while ((indexn != nj) && (iter <= (n - i + 1)))
    {
      nj = ni
      ni = indexn
      indexn = which.min(distmatn[, ni])
      iter = iter + 1         # control for bad trees
    }
    tempindex = c(ni, nj)
    namenew = n+i
    topo[i, 1:2] = name[tempindex]   # 1st & 2nd columns: nodes combined
    topo[i, 3] = namenew    # 3rd column: the number of the new node

    nnodes1 = nnodes[name[ni]]
    nnodes2 = nnodes[name[nj]]
    nnodes[n+i] = nnodes1 + nnodes2
    topo[i, 4] = nnodes[n+i]    # 4th column: the number of the leaves in this new node

    delta1 = mean(distmatn[, ni], na.rm = TRUE)   # delta1 = (sum_k a_{k,i}) / (r-1)
    delta2 = mean(distmatn[, nj], na.rm = TRUE)   # delta2 = (sum_k a_{k,j}) / (r-1)

    dijsum = distmatn[nj, ni] + distmatn[ni, nj]   # sum of aij and aji
    di2node = (delta1 - delta2) / 2 + dijsum / 2   # Pi to new node
    dj2node = (delta2 - delta1) / 2 + dijsum / 2   # Pj to new node

    leng[i, ] = c(di2node, dj2node)

    # update node -> others
    dnew2k = (distmatn[, ni] + distmatn[, nj] - di2node - dj2node) / 2
    distmatn[, nj] = dnew2k
    dk2new = (distmatn[ni, ] + distmatn[nj, ]) / 2
    distmatn[nj, ] = dk2new
    name[nj] = namenew
    distmatn = distmatn[-ni, -ni]
    name = name[-ni]
    #    distmatn = pmax(distmatn, 0)       # make all branch lengthes positive
  }
  list(topo = topo, leng = leng)
}



# treeForm
# generate a tree string in the newick format, e.g. ((A:0.5, B:0.6):0.1, C:0.2);

treeForm = function(topo, leng, name)
{
  leng = round(leng, digits = 5)   # use only the first 5 digits

  n = dim(leng)[1] + 1
  temp = list()
  temp[1:n] = name

  for (i in 1:(n-1))
  {
    temp[n+i] = paste('(', temp[topo[i, 1]], ':', leng[i, 1], ',',
                      temp[topo[i, 2]], ':', leng[i, 2], ')', sep = '')
  }
  res = paste(temp[n+i][[1]], ';', sep = '')
  res
}


# treePlot
# generate a tree plot

treePlot = function(distMat, name)
{
  tem = tree0(distMat)
  res = treeForm(tem$topo, tem$leng, name)
  res = read.tree(text = res)
  plot(res)
}


#####################################################
# Simulation.

# sim2pop
# use true mf
# multiple simulation runs for 2 populations, p0 -> p1, p2

sim2pop = function(p0, nsim, s2q1, s2q2, sigmastart = c(0.5, 0.5),
                   mf = 0.2, zero = zero40, weight = weight40, steptol = 1.e-6)
{
  n = length(p0)

  s2qhat = matrix(NA, nsim, 2)

  for (i in 1:nsim)
  {
    cat(i, '\n')
    p1 = simNormAppr(p0, s2q1)
    p2 = simNormAppr(p0, s2q2)

    p1 = symm(p1)
    p2 = symm(p2)

    res12 = nlm(nllksigmavec, sigmastart, p1 = p1, p2 = p2, zero = zero,
                weight = weight, mf = mf, iterlim = 200, steptol = steptol,
                print.level = 0)

    s2q1hat = res12$est[1]
    s2q2hat = res12$est[2]
    s2qhat[i, ] = c(s2q1hat, s2q2hat)
    cat(s2qhat[i, ], '\n')
  }
  s2qhat
}


# sim2pop1
# estimate mf
# multiple simulation runs for 2 populations, p0 -> p1, p2

sim2pop1 = function(p0, nsim, s2q1, s2q2, sigmastart = c(0.5, 0.5),
                    zero = zero40, weight = weight40, steptol = 1.e-6)
{
  n = length(p0)

  #  quad = gauss.quad(nquad, kind = "legendre")
  #  zero = quad$nodes
  #  weight = quad$weights

#  s2qhat = matrix(NA, nsim, 2)
  s2qhat = matrix(NA, nsim, 3)

  for (i in 1:nsim)
  {
#    cat(i, '\n')
    p1 = simNormAppr(p0, s2q1)
    p2 = simNormAppr(p0, s2q2)

    p1 = symm(p1)
    p2 = symm(p2)

    datat = (p1 + p2)/2
    index.fix0 = datat == 0
    index.fix1 = datat == 1
    m0hat = sum(index.fix0) / length(p1)
    m1hat = sum(index.fix1) / length(p1)
    mfhat = m0hat + m1hat

    res12 = nlm(nllksigmavec, sigmastart, p1 = p1, p2 = p2, zero = zero,
                weight = weight, mf = mfhat, iterlim = 200,
                print.level = 0, steptol = steptol)

    s2q1hat = res12$est[1]
    s2q2hat = res12$est[2]
#    s2qhat[i, ] = c(s2q1hat, s2q2hat)
    # record mfhat
    s2qhat[i, ] = c(s2q1hat, s2q2hat, mfhat)
#    cat(s2qhat[i, ], '\n')
  }
  s2qhat
}


# condCheck
# check the condition that there exists at least one SNP such that
# p1i==0 and p2i!=0, or p1i==1 and p2i!=1 (condition-j,i)
# and there exists at least one SNP such that
# p1j!=0 and p2j==0, or p1j!=1 and p2j==1 (condition-i,j)
# output: a matrix: n*n, the (i,j) element indicats whether
# condition-i,j is satisfied

condCheck = function(data)
{
  n = dim(data)[1]
  res = matrix(NA, n, n)
  for (i in 1:(n-1))
  {
    p1 = data[i, ]
    for (j in (i+1):n)
    {
      p2 = data[j, ]
      tem1 = any((p1==0)&(p2!=0))
      tem2 = any((p1==1)&(p2!=1))
      res[j, i] = any(tem1, tem2)
      tem1 = any((p2==0)&(p1!=0))
      tem2 = any((p2==1)&(p1!=1))
      res[i, j] = any(tem1, tem2)
    }
  }
  res
}


# condCheck2pop
# check the condition that there exists at least one SNP such that
# p1i==0 and p2i!=0, or p1i==1 and p2i!=1 (condition-j,i)
# and there exists at least one SNP such that
# p1j!=0 and p2j==0, or p1j!=1 and p2j==1 (condition-i,j)
# output: TRUE or FALSE

condCheck2pop = function(p1, p2)
{
  tem1 = any((p1==0)&(p2!=0))
  tem2 = any((p1==1)&(p2!=1))
  res1 = any(tem1, tem2)
  tem1 = any((p2==0)&(p1!=0))
  tem2 = any((p2==1)&(p1!=1))
  res2 = any(tem1, tem2)
  res1 & res2
}


# condNum
# count the condition numbers that there exists at least one SNP such that
# p1i==0 and p2i!=0, or p1i==1 and p2i!=1 (condition-j,i)
# and there exists at least one SNP such that
# p1j!=0 and p2j==0, or p1j!=1 and p2j==1 (condition-i,j)
# output: a matrix: n*n, the (i,j) element indicats whether
# condition-i,j is satisfied

condNum = function(data)
{
  n = dim(data)[1]
  res = matrix(NA, n, n)
  for (i in 1:(n-1))
  {
    p1 = data[i, ]
    for (j in (i+1):n)
    {
      p2 = data[j, ]
      tem1 = sum((p1==0)&(p2!=0))
      tem2 = sum((p1==1)&(p2!=1))
      res[j, i] = tem1 + tem2
      tem1 = sum((p2==0)&(p1!=0))
      tem2 = sum((p2==1)&(p1!=1))
      res[i, j] = tem1 + tem2
    }
  }
  res
}


# simTree
# simulate a data set of a given tree topology, using truncated Normal model
# input: tree is a phylogenetic tree

simTree = function(tree, nsnp = 1000)
{
  # collect information from tree
  nnode = tree$Nnode
  nleaf = nnode+1
  label = tree$tip.label
  edge = tree$edge
  leng = tree$edge.length

  # initialize p0, the freqs at the root
  #  p0 = runif(nsnp)

  p0 = c(runif(nsnp * 0.8), rep(c(0, 1), each = nsnp * 0.1))

  freq = matrix(NA, 2*nleaf-1, nsnp)   # all freqs, including the internal nodes and leaves
  freq[nleaf+1, ] = p0    # initialize p0, which is by nature the (nleaf+1) row
  for (i in 1:(2*nleaf-2))
  {
    index1 = edge[i, 1]
    index2 = edge[i, 2]
    p0 = freq[index1, ]
    p1 = simNormAppr(p0, leng[i])    # simulate using the trunNorm model
    freq[index2, ] = p1
  }
  freq = freq[1:nleaf, ]
  dimnames(freq) = list(label)
  freq
}


# treeSimEst

# simulate n data set of a given tree topology, using truncated Normal model
# and return estimated trees for each
# input: tree is a phylogenetic tree

treeSimEst = function(treeTrue, nsim = 10, nsnp = 1000, zero = zero40,
                      weight = weight40, steptol = 1.e-6, sigmastart = c(0.2, 0.2))
{
  nnode = length(treeTrue$Nnode)

  # result using our method
  res = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)

  # result using
  resnj = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)
  for (i in 1:nsim)
  {
    cat("\n", "          nsim = ", i, "           ", "\n")
    treeSimFreq = simTree(treeTrue, nsnp)
    treeMat = distm1(treeSimFreq, sigmastart, zero, weight, steptol)
    treeSimEst = tree0(treeMat)
    treeSimEstForm = treeForm(treeSimEst$topo, treeSimEst$leng,
                              rownames(treeSimFreq))
    treeSimEst = read.tree(text = treeSimEstForm)
    res[[i]] = treeSimEst
    resnj[[i]] = nj(dist.prop(as.data.frame(treeSimFreq), 4))
  }
  list(res = res, resnj = resnj)
}


# treeSimEstNJ

# only use symmetric dissimilarity matrices, but different distance measures
# simulate n data set of a given tree topology, using truncated Normal model
# and return estimated trees for each
# input: tree is a phylogenetic tree

treeSimEstNJ = function(treeTrue, nsim = 10, nsnp = 1000, sigmastart = c(0.2, 0.2),
                        zero = zero40, weight = weight40)
{
  nnode = length(treeTrue$Nnode)

  # result using dist 1 in dist.prop
  resnj1 = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)
  # result using dist 1 in dist.prop
  resnj2 = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)
  # result using dist 1 in dist.prop
  resnj3 = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)
  # result using dist 1 in dist.prop
  resnj4 = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)
  # result using dist 1 in dist.prop
  resnj5 = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)

  for (i in 1:nsim)
  {
    cat("\n", "          nsim = ", i, "           ", "\n")
    treeSimFreq = simTree(treeTrue)
    resnj1[[i]] = nj(dist.prop(as.data.frame(treeSimFreq), 1))
    resnj2[[i]] = nj(dist.prop(as.data.frame(treeSimFreq), 2))
    resnj3[[i]] = nj(dist.prop(as.data.frame(treeSimFreq), 3))
    resnj4[[i]] = nj(dist.prop(as.data.frame(treeSimFreq), 4))
    resnj5[[i]] = nj(dist.prop(as.data.frame(treeSimFreq), 5))
  }
  list(resnj1 = resnj1, resnj2 = resnj2, resnj3 = resnj3, resnj4 = resnj4,
       resnj5 = resnj5)
}


# scaleTree
# scale a tree so that the total length sum up to 1 (by default)
# input: tree is a phylogenetic tree

scaleTree = function(tree, tlength = 1)
{
  # collect information from tree
  leng = tree$edge.length
  leng = leng / sum(leng) * tlength
  tree$edge.length = leng
  tree
}


# consensus tree with branch lengths from those with correct topology

conTreeEdgeLen = function(trees)
{
  ntree = length(trees)
  resTree = trees[[1]]
  nedge = length(resTree$edge.length)
  temp = matrix(NA, ntree, nedge)
  for (i in 1:ntree)
  {
    temp[i, ] = trees[[i]]$edge.length
  }
  resTree$edge.length = apply(temp, 2, mean)
  resTree
}


# treeDictRoot
# dictionary form of tree representation
# for one tree

treeDictRoot = function(tree)
{
  nnode = tree$Nnode
  nedge = length(tree$edge.length)
  keynames = rep(NA, nedge+1)
  ntip = length(tree$tip.label)
  keynames[1:ntip] = tree$tip.label
  subbranch = list(NA, nedge+1)
  sublength = list(NA, nedge+1)
  known = 1:ntip
  unknown = (ntip+1):(nedge+1)
  lenU = length(unknown)
  # nnode = length(unknown)
  for (i in 1:lenU)
  {
    for (j in 1:(lenU-i+1))
    {
      temnode = unknown[j]
      temtipsindex = tree$edge[, 1] == temnode
      temtips = tree$edge[temtipsindex, 2]
      if (!all(temtips %in% known))
      { next } else
      {
        temname = keynames[temtips]
        temname1 = unlist(strsplit(temname[1], '_'))
        temname2 = unlist(strsplit(temname[2], '_'))
        temnamevec = c(temname1, temname2)
        ordername = order(temname)
        temnamenew = sort(temnamevec)
        keynames[temnode] = paste(temnamenew, collapse = '_')
        temedgelength = tree$edge.length[temtipsindex]
        sublength[[temnode]] = temedgelength[ordername]
        subbranch[[temnode]] = temname[ordername]
        known = c(known, temnode)
        unknown = unknown[unknown != temnode]
        break
      }
    }
  }
  hashkeyname = keynames[(ntip+1):(nedge+1)]
  hashlength = sublength[(ntip+1):(nedge+1)]
  hashbranch = subbranch[(ntip+1):(nedge+1)]
  reslength = hash(keys = hashkeyname, values = hashlength)
  resbranch = hash(keys = hashkeyname, values = hashbranch)
  subjointbranch = rep(NA, length(hashbranch))
  for (i in 1:length(hashbranch))
  {
    subjointbranch[i] = paste(hashbranch[[i]], collapse = " & ")
  }
  allbranch = paste(hashkeyname, subjointbranch, sep = " -> ")
  resboth = hash(keys = allbranch, values = hashlength)
  list(dictBranch = resbranch, dictLength = reslength, dictBoth = resboth)
}


# treeDictUnroot
# dictionary form of tree representation
# for one tree

treeDictUnroot = function(tree)
{
  nnode = tree$Nnode
  nedge = length(tree$edge.length)
  keynames = rep(NA, nedge+1)
  ntip = length(tree$tip.label)
  keynames[1:ntip] = tree$tip.label
  subbranch = list(NA, nedge+1)
  sublength = list(NA, nedge+1)
  known = 1:ntip
  unknown = (ntip+1):(nedge+1)
  lenU = length(unknown)
  for (i in 1:(lenU-1))
  {
    for (j in 1:(lenU-i+1))
    {
      temnode = unknown[j]
      temtipsindex = tree$edge[, 1] == temnode
      temtips = tree$edge[temtipsindex, 2]
      if (!all(temtips %in% known))
      { next } else
      {
        temname = keynames[temtips]
        temname1 = unlist(strsplit(temname[1], '_'))
        temname2 = unlist(strsplit(temname[2], '_'))
        temnamevec = c(temname1, temname2)
        ordername = order(temname)
        temnamenew = sort(temnamevec)
        keynames[temnode] = paste(temnamenew, collapse = '_')
        temedgelength = tree$edge.length[temtipsindex]
        sublength[[temnode]] = temedgelength[ordername]
        subbranch[[temnode]] = temname[ordername]
        known = c(known, temnode)
        unknown = unknown[unknown != temnode]
        break
      }
    }
  }
  temnode = unknown[1]
  temtipsindex = tree$edge[, 1] == temnode
  temtips = tree$edge[temtipsindex, 2]

  temname = keynames[temtips]
  temname1 = unlist(strsplit(temname[1], '_'))
  temname2 = unlist(strsplit(temname[2], '_'))
  temname3 = unlist(strsplit(temname[3], '_'))
  temnamevec = c(temname1, temname2, temname3)
  ordername = order(temname)
  temnamenew = sort(temnamevec)
  keynames[temnode] = paste(temnamenew, collapse = '_')
  temedgelength = tree$edge.length[temtipsindex]
  sublength[[temnode]] = temedgelength[ordername]
  subbranch[[temnode]] = temname[ordername]

  hashkeyname = keynames[(ntip+1):(nedge+1)]
  hashlength = sublength[(ntip+1):(nedge+1)]
  hashbranch = subbranch[(ntip+1):(nedge+1)]
  reslength = hash(keys = hashkeyname, values = hashlength)
  resbranch = hash(keys = hashkeyname, values = hashbranch)
  subjointbranch = rep(NA, length(hashbranch))
  for (i in 1:length(hashbranch))
  {
    subjointbranch[i] = paste(hashbranch[[i]], collapse = " & ")
  }
  allbranch = paste(hashkeyname, subjointbranch, sep = " -> ")
  resboth = hash(keys = allbranch, values = hashlength)
  list(dictBranch = resbranch, dictLength = reslength, dictBoth = resboth)
}


# treeDict
# for one tree
# works for both rooted trees and unrooted trees

treeDict = function(tree)
{
  if (is.rooted(tree))
  {
    res = treeDictRoot(tree)
  } else
  {
    res = treeDictUnroot(tree)
  }
  res
}


# treesDict
# for n trees
# dictionary form of tree representation

treesDict = function(trees)
{
  ntree = length(trees)
  reslength = list()
  resbranch = list()
  resboth = list()
  for (i in 1:ntree)
  {
    temp = treeDict(trees[[i]])
    resbranch[[i]] = temp$dictBranch
    reslength[[i]] = temp$dictLength
    resboth[[i]] = temp$dictBoth
  }
  list(dictBranch = resbranch, dictLength = reslength, dictBoth = resboth)
}


# treeCon
# consensus tree of all nodes
# using dictionary form of trees
# treesDictFrom : dictBranch from treesDict
# output: similar to prop.part in package ape

treeCon = function(treesDictFormBranch)
{
  trees = treesDictFormBranch
  ntree = length(trees)

  # initialize
  branchall = keys(trees[[1]])
  for (i in 2:ntree)
  {
    branchall = c(branchall, keys(trees[[i]]))
  }

  uBranch = unique(branchall)
  nuBranch = length(uBranch)
  nBranch = rep(NA, nuBranch)
  for (i in 1:nuBranch)
  {
    branch = uBranch[i]
    branchindex = branchall == branch
    nBranch[i] = sum(branchindex)
  }
  nBranchCon = hash(keys = uBranch, values = nBranch)
  nBranchCon
}


# treeConOut
# output, a consensus tree with branch lengths

treeConOut = function(treeConExactRes)
{
  treeConNum = treeConExactRes$treeConNum
  treeConLen = treeConExactRes$treeConLen
  valuesNum = values(treeConNum)
  temp = keys(treeConNum)[valuesNum < 0.5*max(valuesNum)]
  if (length(temp) > 0)
  { del(temp, treeConLen) }
  treeConLen
}


# treeConRes

treeConRes = function(trees, fun='mean')
{
  treesDictBoth = treesDict(trees)$dictBoth
  treeConRes = treeConExact(treesDictBoth, fun)
  treeConRes = treeConOut(treeConRes)
  treeConRes = treeDict2Newick(treeConRes)
  treeConRes
}


# distTree
# distance between tree tree and simulated trees
# output: ntree * 3 matrix

distTree = function(treeTrue, treeSim)
{
  ntree = length(treeSim)
  distres = matrix(NA, ntree, 5)
  for (i in 1:ntree)
  {
    treei = treeSim[[i]]
    treeiScaled = scaleTree(treei, sum(treeTrue$edge.length))
    distres[i, 1] = dist.topo(treeTrue, treei)
    distres[i, 2] = dist.topo(treeTrue, treei, "score")
    distres[i, 3] = dist.topo(treeTrue, treeiScaled, "score")
    distres[i, 4] = dist.multiPhylo(c(treeTrue, treei))
    distres[i, 5] = dist.multiPhylo(c(treeTrue, treeiScaled))
  }
  #  rownames(distres) = c('PH', 'KF', 'KF-scaled', 'geodesic', 'geodesic-scaled')
  distres
}


# treeConExact
# exact consensus pool of all nodes
# i.e., correct nodes with correct sub nodes at next level
# input: using dictionary form of trees
# output: number of correct cherries (or internal cherries)
#         median branch length of those cherries

treeConExact = function(treesDictFormBoth, fun='mean')
{
  trees = treesDictFormBoth
  ntree = length(trees)

  # initialize
  branchall = keys(trees[[1]])
  steplen = length(branchall)
  lengthall = list(NA, ntree * steplen)
  lengthall[1:steplen] = values(trees[[1]], simplify = F)
  for (i in 2:ntree)
  {
    branchall = c(branchall, keys(trees[[i]]))
    lengthall[((i-1)*steplen+1):(i*steplen)] = values(trees[[i]], simplify = F)
  }

  uBranch = unique(branchall)
  nuBranch = length(uBranch)
  nBranch = rep(NA, nuBranch)
  lBranch = list(NA, nuBranch)
  for (i in 1:nuBranch)
  {
    branch = uBranch[i]
    branchindex = branchall == branch
    nBranch[i] = sum(branchindex)
    temlen = matrix(unlist(lengthall[branchindex]), nBranch[i],
                    length(unlist(lengthall[branchindex])) / nBranch[i],
                    byrow = TRUE)
    if (fun=='var')
    {
      lBranch[[i]] = apply(temlen, 2, var)
    }
    else
    {
      lBranch[[i]] = apply(temlen, 2, mean)      
    }
  }
  nBranchCon = hash(keys = uBranch, values = nBranch)
  lBranchCon = hash(keys = uBranch, values = lBranch)
  list(treeConNum = nBranchCon, treeConLen = lBranchCon)
}


# treeDict2Newick
# tree changing from dictionary form to Newick form
# input: treeDictBoth format

treeDict2Newick = function(treeDictBoth)
{
  tree = treeDictBoth
  nnode = length(tree)
  branch = keys(tree)
  leng = values(tree, simplify = F)

  parent = rep(NA, nnode)
  nparent = rep(NA, nnode)
  child = list(NA, nnode)
  nchild = matrix(NA, nnode)
  for (i in 1:nnode)
  {
    temp = unlist(strsplit(branch[i], " -> "))
    parent[i] = temp[1]
    nparent[i] = length(unlist((strsplit(temp[1], '_'))))
    child[[i]] = unlist(strsplit(temp[2], " & "))
  }
  orderp = order(nparent, decreasing = T)
  parent = parent[orderp]
  nparent = nparent[orderp]
  child = child[orderp]
  leng = leng[orderp]

  childlength = rep(NA, length(parent))
  if (length(child[[1]]) == 3)
  {
    child1length1 = paste(child[[1]][1], leng[[1]][1], sep = ':')
    child1length2 = paste(child[[1]][2], leng[[1]][2], sep = ':')
    child1length3 = paste(child[[1]][3], leng[[1]][3], sep = ':')
    child1length = paste(child1length1, child1length2, child1length3, sep = ',')
    childlength[1] = paste('(', child1length, ')', sep = '')
  } else
  {
    child1length1 = paste(child[[1]][1], leng[[1]][1], sep = ':')
    child1length2 = paste(child[[1]][2], leng[[1]][2], sep = ':')
    child1length = paste(child1length1, child1length2, sep = ',')
    childlength[1] = paste('(', child1length, ')', sep = '')
  }
  for (i in 2:length(parent))
  {
    childlength1 = paste(child[[i]][1], leng[[i]][1], sep = ':')
    childlength2 = paste(child[[i]][2], leng[[i]][2], sep = ':')
    childlengthi = paste(childlength1, childlength2, sep = ',')
    childlength[i] = paste('(', childlengthi, ')', sep = '')
  }

  res = paste(parent[1], ';', sep = '')
  for (i in 1:length(parent))
  {
    temparent = parent[i]
    res = sub(temparent, childlength[i], res)
  }
  read.tree(text = res)
}


# simTreeWF
# simulate data based ona a tree and Wright-Fisher model

simTreeWF = function(tree, nsnp = 1000, npopsize=1000)
{
  # collect information from tree
  nnode = tree$Nnode
  nleaf = nnode+1
  label = tree$tip.label
  edge = tree$edge
  leng = tree$edge.length

  # initialize p0, the freqs at the root
  #  p0 = runif(nsnp)

  p0 = c(runif(nsnp * 0.8), rep(c(0, 1), each = nsnp * 0.1))

  freq = matrix(NA, 2*nleaf-1, nsnp)   # all freqs, including the internal nodes and leaves
  freq[nleaf+1, ] = p0    # initialize p0, which is by nature the (nleaf+1) row
  for (i in 1:(2*nleaf-2))
  {
    index1 = edge[i, 1]
    index2 = edge[i, 2]
    p0 = freq[index1, ]
    p1 = simWF(p0, leng[i], N=npopsize)    # simulate using the WF model
    freq[index2, ] = p1
  }
  freq = freq[1:nleaf, ]
  dimnames(freq) = list(label)
  freq
}

treeSimEstNJ = function(treeTrue, nsim = 10, nsnp = 1000, zero = zero40,
                        weight = weight40, steptol = 1.e-6, sigmastart = c(0.2, 0.2))
{
  nnode = length(treeTrue$Nnode)

  # result using our method
  res = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)

  # result using
  resnj = rmtree(nsim, nnode + 1, rooted = TRUE, tip.label = NULL)
  for (i in 1:nsim)
  {
    cat("\n", "          nsim = ", i, "           ", "\n")
    treeSimFreq = simTreeWF(treeTrue, nsnp)
    treeMat = distm1(treeSimFreq, sigmastart, zero, weight, steptol)
    treeSimEst = tree0(treeMat)
    treeSimEstForm = treeForm(treeSimEst$topo, treeSimEst$leng,
                              rownames(treeSimFreq))
    treeSimEst = read.tree(text = treeSimEstForm)
    res[[i]] = treeSimEst
    resnj[[i]] = nj(dist.prop(as.data.frame(treeSimFreq), 4))
  }
  list(res = res, resnj = resnj)
}


# simWF
# simulation based on the Wright-Fisher model
# p: initial population frequence
# nsim: number of simulation runs
# t: int, generations
# N: population size, fixed for the process

simWF = function(p, t, N = 1000)
{
  nsim = length(p)
  res = p
  pnew = res
  for (i in 1:t)
  {
    xnew = rbinom(nsim, N, pnew)
    pnew = xnew / N
    #  	res = rbind(res, pnew)    # show the history of p
  }
  pnew
}

