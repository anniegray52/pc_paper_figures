library("transport")
library("dimension")
library("igraph")
library("TDA")
library(Matrix)
library(dimRed)
library("cccd")

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

## ladle method commented out for speed

## ================= initial settings

split = .5
minpn = 500
sdeps=1
qepsilon = .05
qmaxscale = .2
bsub = min(minpn, 200)

## ==================== functions 


mdist = function(Y, epsilon){
    DistY = as.matrix(dist(Y))
    newdist=DistY
    newdist[DistY > epsilon] = 0
    epsilongraph = graph_from_adjacency_matrix(newdist, weighted=TRUE,mode="min")
    ManifoldDist = distances(epsilongraph); ManifoldDist[ManifoldDist==Inf]=-1
    pmax(ManifoldDist, DistY)
}

scale = function(x, benchmark){((x-min(x))/(max(x)-min(x)))*(max(benchmark)-min(benchmark))  + min(benchmark)}

standardise <- function(X){
    Xc = X - colMeans(X)
    Xc/(sqrt(mean(rowSums(Xc^2))))
}

rotation = function(theta){
    matrix(c(cos(theta),sin(theta),-sin(theta),cos(theta)),ncol =2, byrow= TRUE)
}

make_circle = function(n,x = 0,y = 0, r = 1){
    s = runif(n, 0, 2*pi)
    cbind(r*cos(s) + x, r*sin(s) + y)
}

create_Z = function(n= 1000, a = -1, b = 1, c = -1, d = 1, q = .5){
    r = q * tan(pi/8)
    theta = 2 * atan(1/tan(pi/8))
    e = q * sin(pi/4)
    xstart = a - sign(a)* e
    xend = b - sign(b)* e
    ystart = ((d-c)/(b-c)) * xstart
    yend = ((d-c)/(b-c)) * xend
    sll1 = (b - sign(b)*q) - (a - sign(a)*(q-r))
    sll3 = sqrt((yend-ystart)**2 + (xend - xstart)**2)
    cl1 = r*theta
    ls = c(sll1,sll1,sll3,cl1,cl1)
    total_l = sum(ls)
    sizes <- c(rmultinom(n = 1, size = n, prob = ls/total_l))
    
    sl1 = cbind(runif(sizes[1], a - sign(a)*(q-r) , b - sign(b)*q), rep(d,sizes[1]))
    sl2 = cbind(runif(sizes[2], a - sign(a)*q , b - sign(b)*(q-r)), rep(c,sizes[2]))
    
    x = runif(sizes[3], xstart, xend)
    y = ((d-c)/(b-c)) * x 
    sl3 = cbind(x,y)
    sls = rbind(sl1,sl2,sl3)
    
    c1 = runif(sizes[4],0, theta)
    c2 = runif(sizes[5],0, theta)
    
    cs1 = cbind(cos(c1),sin(c1)) %*% rotation(-pi/4)
    cs2 = cbind(cos(c2),sin(c2)) %*% rotation(3*pi/4)
    cs = rbind(cbind(cs1[,1]*r + b-q, cs1[,2]*r + d-r),cbind(cs2[,1]*r + a+q, cs2[,2]*r + c+r))
    
    rbind(sls,cs)
}


computeK <- function(Z) {
    K = matrix(0, nrow(Z), nrow(Z))
    for (i in 1:nrow(Z)) {
        for (j in 1:nrow(Z)) {
            K[i, j] = kernel(Z[i, ], Z[j, ])
        }
    }
    K
}


## ======================== run


experiment_results = list()
for (np in 1:2) {
    for (kr in 1:4) {
        print("Here")
        {
            # np = 1
            # kr = 4
            
            if (np == 1) {
                n = minpn
                p = 2 * n
            } else {
                p = minpn
                n = 2 * p
            }
            if (kr == 1) {
                ##K points
                bd = 0
                nc = 6
                Z = (1:nc)[sample(nc, n, replace = TRUE)]
                dim(Z) = c(n, 1)
                kernel = function(x, y)
                    (exp(-1 / 2 * sum((
                        x - y
                    ) ^ 2)))
            }
            if (kr == 2) {
                ##circles, polynomial kernel,
                bd = 1
                r = 2
                kernel = function(x, y)
                    (t(x) %*% y + 1) ^ r
                # theta = runif(n,-pi, pi)
                # Z = cbind(cos(theta), sin(theta))
                sizes = c(rmultinom(n = 1, size = n, prob = rep(1/4, 4)))
                Z = rbind(make_circle(sizes[1],0,1),make_circle(sizes[2],0,-1),make_circle(sizes[3],1,0),make_circle(sizes[4],-1,0))
            }
            if (kr == 3) {
                ##Z, cosine kernel
                bd = 0
                kernel = function(x, y) {
                    (cos(x[1] - y[1]) + cos(x[2] - y[2]) + 2)
                }
                # Z = cbind(runif(n,-pi + 0.25, pi - 0.25),
                #           runif(n,-pi + 0.25, pi - 0.25))
                Z = create_Z(n)
            }
            if (kr == 4) {
                bd = 1
                ##circle, RBF kernel,
                kernel = function(x, y)
                    (exp(-1 / 2 * sum((
                        x - y
                    ) ^ 2)))
                theta = runif(n,-pi, pi)
                Z = cbind(cos(theta), sin(theta))
            }
            
            K = computeK(Z)
            E = eigen(K)
            ####checking it's positive definite
            ##round(E$values, digits=5)
            # rank = sum(round(E$values, digits = 10) > 0)
            rank = rankMatrix(K)[1]
            print(rank)
            L = E$vectors[, 1:rank] %*% diag(sqrt(E$values[1:rank]))
            N = matrix(rnorm(rank * p), nrow = rank, ncol = p)
            X = L %*% N
            ####X has the same rank as K
            ##rankX = sum(round(svd(X)$d, digits=5)>0)
            
            DistZ = as.matrix(dist(Z))
            epsilon = quantile(DistZ, qepsilon)
            maxscale = 1
            TrueM = mdist(Z, epsilon)
            ##Zs = standardise(Z);
            Zs=Z
            tdiagram = (ripsDiag(
                Zs[1:bsub, ],
                maxdimension = 1,
                maxscale = maxscale
            ))$diagram
            
            Y = X + matrix(rnorm(n * p, mean = 0, sd = sdeps),
                           nrow = n,
                           ncol = p)
            ###dimension selection
            train = round(n * split)
            
            rmax = min(train, 50) ##was 100
            YPCA = 1/sqrt(p)*(prcomp(Y, center = FALSE, rank = rmax))$x
             
            Ytrain = Y[1:train, ]
            S = svd(Ytrain)
            ##plot(S$d)
            ws = c()
            projdistances = c()
            geodesic_error = c()
            bottleneck_error = c()
            tda_error = c()
            for (rtry in 1:rmax) {
                P = S$v[, 1:rtry, drop = FALSE] %*% t(S$v[, 1:rtry, drop = FALSE])
                Yproj = Ytrain %*% t(P)
                w = transport::wasserstein(pp(Yproj), pp(Y[(train + 1):nrow(Y), ]))
                ##print(w)
                ws = c(ws, w)
                projdistance = sum((Y[(train + 1):nrow(Y), ] - Y[(train + 1):nrow(Y), ] %*% t(P)) ^
                                       2)
                projdistances = c(projdistances, projdistance)
                
                geodesic_error = c(geodesic_error, mean((
                    TrueM - mdist(YPCA[, 1:rtry], epsilon)
                ) ^ 2))
                ##YPCAs = standardise(YPCA[, 1:rtry, drop = FALSE])
                YPCAs = YPCA[, 1:rtry, drop = FALSE]
                sdiagram = (ripsDiag(
                    YPCAs[1:bsub, ],
                    maxdimension = 1,
                    maxscale = maxscale
                ))$diagram
                bottleneck_error = c(bottleneck_error,
                                     bottleneck(tdiagram, sdiagram, dimension = bd))
            }
            # res=dimension(Y, method="ladle")
            experiment_results[[length(experiment_results) + 1]] = list(
                Z = Z,
                YPCA = YPCA,
                rank = rank,
                X = X,
                Y = Y,
                ws = ws,
                projdistances = projdistances,
                geodesic_error = geodesic_error,
                bottleneck_error = bottleneck_error
                # ,ladle=res
            )
        }
        
    }
}

## save.image("dimension_selection_comparison_minpn500.RData")
## load("dimension_selection_comparison_minpn500.RData")

## ##plot(prcomp(Y,center=FALSE)$x[,1:2])
## main=c("circle, polynomial kernel", "circle, RBF", "{1, ..., 10}, RBF")
##dev.new(width=10, height=15)

e <-function(){rnorm(1,0,.2)}
pdf("dimension_selection_experiments.pdf", width=10, height=13)
par(mfrow=c(5,4))
par(mar=c(4.5, 4.5, 2.5, 1.5))
for (kr in 1:4){
    results = experiment_results[[kr]];
    if (kr > 1){
        plot(results$Z, xlab="", ylab="", pch=16, cex=.5)
    }
    else {        
        plot(cbind(rep(1:2, 3), rep(1:3, 2)), xlim = c(0,3), ylim=c(0.5,3.5), xlab="", ylab="", pch=16, cex=.5, axes="F")
    }
    if (kr==1){mtext("a", side=3, adj=0, line=0.5)}
    mtext(paste0(kr,""), side=3, line=0.5)
}

##PC1, PC2
for (kr in 1:4){
    results = experiment_results[[kr]];
    plot(results$YPCA[,1:2], xlab="PC1", ylab="PC2", pch=16, cex=.5)
    if (kr==1){mtext("b", side=3, adj=0, line=0.5)}
}


##Dimension selection
for (kr in 1:4){
    results = experiment_results[[kr]];
    plot(log(results$ws), xlim = c(1,30), ylab="", xlab="rank", yaxt="n");
    if (kr==1){mtext("c", side=3, adj=0, line=0.5)}
    title(ylab="Wasserstein error", line=1)
        abline(v=which.min(results$ws), col="red");
        abline(v=dim_select((svd(results$Y))$d)+e(), col="orange", lwd=1)
        # abline(v=which.min((results$ladle)$gn)-1+e(), col="blue", lwd=1)
    if (kr == 1) {
                                        #legend("topright", legend=c("true rank", "Wasserstein method", "Elbow method", "Ladle method", "Wasserstein distance"), col=c("black", "red", "orange", "blue", "black"), lty=c(2,1,1,1,NA), pch= c(NA,NA,NA,NA,1))}
        legend("topright", legend=c("true rank", "Wasserstein", "Elbow", "Ladle"), col=c("black", "red", "orange", "blue"), lty=c(2,1,1,1))}
        if (kr != 4) {
            abline(v=results$rank, col="black", lwd=2, lty=2);
        }
        ## if (kr >=3){
        ##     lines(scale(log(results$geodesic_error), log(results$ws)), col="blue")
        ## }

        ## lines(scale(log(results$bottleneck_error), log(results$ws)), col="purple")
}

##Geodesic error
for (kr in 1:4){
    if (kr <=2){plot.new()}
    else {
    results = experiment_results[[kr]];
    plot(log(results$ws), xlim = c(1,30), ylab="", xlab="rank", yaxt="n");
    title(ylab="Error", line=1)
    if (kr == 3) {
        legend("bottomright", legend=c("Wasserstein", "Geodesic", "Bottleneck"), col=c("black", "blue", "red"), lty=c(NA,1,1), pch= c(1,NA,NA))}
        ## if (kr != 4) {
        ##     abline(v=results$rank + e(), col="black", lwd=2, lty=2);
        ## }
    lines(scale(log(results$geodesic_error), log(results$ws)), col="blue")
    lines(scale(log(results$bottleneck_error), log(results$ws)), col="red")
    }
    if (kr==1){mtext("d", side=3, adj=0, line=0.5)}
}

##bsub=100
bsub=500
pds = list()
# load("persistent_diagrams_dim_selection.Rout")
for (kr in 1:4){
    results = experiment_results[[kr]];
    rank = which.min(results$ws)
    YPCA = results$YPCA[, 1:rank, drop = FALSE]
    sdiagram = (ripsDiag(
        YPCA[1:bsub, ],
        maxdimension = 1,
        maxscale = 4))$diagram
    pds[[length(pds)+1]] = sdiagram
    ## plot(pds[[kr]], cex=.5)
    pd=pds[[kr]]
    plot(pd[,2:3], col=pd[,1]+1, pch=16, xlim=range(c(pd[,2:3])), ylim=range(c(pd[,2:3])), cex=.5)
abline(a = 0, b=1, lty=1)
thr=0.2;
abline(a = thr, b=1, lty=2)

pd0=pd[pd[,1]==0,2:3]; beta0 = sum((pd0[,2]-pd0[,1])>=thr)
pd1=pd[pd[,1]==1,2:3]; beta1 = sum((pd1[,2]-pd1[,1])>=thr)
##mtext(bquote(hat(beta)[0] == .(beta0)),side=3, line=-2, cex=.8)
text(x=1.5, y = 3.5, bquote(hat(beta)[0] == .(beta0)), cex=1, adj=0)
text(x=1.5, y = 3, bquote(hat(beta)[1] == .(beta1)), cex=1,adj=0, col="red")
legend("bottomright", legend=c("Dim 0", "Dim 1"), pch=c(16,16), col=c("black", "red"))
    if (kr==1){mtext("e", side=3, adj=0, line=0.5)}
}

## expression(paste(hat[beta]_)), side=3, line=-1.5, cex=.8)
##save(pds, file="persistent_diagrams_dim_selection.Rout")

dev.off()






###Supplementary figures
##rank and variance
n = 500; p = 1000
kernel = function(x, y) (exp(-1/2 * sum((x - y)^2)))
theta = runif(n,-pi, pi)
Z = cbind(cos(theta), sin(theta)) 

K = computeK(Z)
E = eigen(K)
####checking it's positive definite
##round(E$values, digits=5)
                                        # rank = sum(round(E$values, digits = 10) > 0)
rank = rankMatrix(K)[1]
print(rank)
L = E$vectors[, 1:rank] %*% diag(sqrt(E$values[1:rank]))
N = matrix(rnorm(rank * p), nrow = rank, ncol = p)
X = L %*% N
####X has the same rank as K
##rankX = sum(round(svd(X)$d, digits=5)>0)

sdepsvar = c(0, .5, 1, 2)
wsl=list()
for (sdepsv in sdepsvar){
    Y = X + matrix(rnorm(n * p, mean = 0, sd = sdepsv),
               nrow = n,
               ncol = p)

train = round(n * split)
rmax = min(train, 50) ##was 100
Ytrain = Y[1:train, ]
S = svd(Ytrain)
ws = c()

for (rtry in 1:rmax) {
    P = S$v[, 1:rtry, drop = FALSE] %*% t(S$v[, 1:rtry, drop = FALSE])
    Yproj = Ytrain %*% t(P)
    w = transport::wasserstein(pp(Yproj), pp(Y[(train + 1):nrow(Y), ]))
    ws = c(ws, w)    
}
    wsl[[length(wsl)+1]] = ws
}

pdf(file="Wasserstein_varying_sigma.pdf", width=5, height=5)
plot(log(wsl[[1]]), xlim = c(1,30), xlab="rank", yaxt="n", ylab="Wasserstein error", type="l");
lines(scale(log(wsl[[2]]),log(wsl[[1]])), xlim = c(1,30), col="red");
lines(scale(log(wsl[[3]]),log(wsl[[1]])), xlim = c(1,30), col="blue");
lines(scale(log(wsl[[4]]),log(wsl[[1]])), xlim = c(1,30), col="orange");
legend("topright", legend=c(expression(paste(sigma, " = 0")), expression(paste(sigma, " = 0.5")), expression(paste(sigma, " = 1")), expression(paste(sigma, " = 2"))), lty=1, col=c("black", "red", "blue", "orange"))
dev.off()




##indistinguishability
pdf("indistinguishability.pdf", width=10, height=2.7)
par(mfrow=c(1,4))
par(mar=c(4.5, 4.5, 2.5, 1.5))
chosenD = rbind(c(3,2), c(3,2), c(3,2), c(3,2))
##PC1, PC2
for (kr in 1:4){
    results = experiment_results[[kr]];
    ##estimated_rank=5    
    ##plot(results$YPCA[,(estimated_rank-1):estimated_rank], xlab=paste0("PC", estimated_rank-1), ylab=paste0("PC", estimated_rank), pch=16, cex=.5)
    plot(results$Y[,1:2], xlab="Y1", ylab="Y2", pch=16, cex=.5)
    ##mtext(paste0(kr,""), side=3, line=0.5)
    ## if (kr==1){mtext("b", side=3, adj=0, line=0.5)}
}
dev.off()



##geometry in higher dimensions
pdf("PCA_other_dimensions.pdf", width=10, height=2.7)
par(mfrow=c(1,4))
par(mar=c(4.5, 4.5, 2.5, 1.5))
chosenD = rbind(c(3,2), c(3,2), c(3,2), c(3,2))
##PC1, PC2
for (kr in 1:4){
    results = experiment_results[[kr]];
    estimated_rank = which.min(results$ws)
    ##estimated_rank=5    
    ##plot(results$YPCA[,(estimated_rank-1):estimated_rank], xlab=paste0("PC", estimated_rank-1), ylab=paste0("PC", estimated_rank), pch=16, cex=.5)
    plot(results$YPCA[,chosenD[kr,]], xlab=paste0("PC", chosenD[kr,1]), ylab=paste0("PC", chosenD[kr,2]), pch=16, cex=.5)
    ##mtext(paste0(kr,""), side=3, line=0.5)
    ## if (kr==1){mtext("b", side=3, adj=0, line=0.5)}
}
dev.off()

par(mfrow=c(4,4))
par(mar=c(4.5, 4.5, 2.5, 1.5))
results = experiment_results[[4]];
for (dim1 in 1:4){
    for (dim2 in 1:4){
        if (dim1 == dim2) plot.new()
        else     plot(results$YPCA[,c(dim1,dim2)], xlab=paste0("PC", dim1), ylab=paste0("PC", dim2), pch=16, cex=.5)
    }        
}
