library(data.table)
library(GMMAT)
library(doMC)
library(foreach)

sink(snakemake@log[[1]])

all.Rdata = snakemake@input[['data']]

all.cohort.scores = list()
all.cohort.covms = list()
cohort.names = NULL
for (Rdata in all.Rdata) {
    load(Rdata)
    cohort_name = sub(".*/([^/]+)/chr[0-9XYM]+/.*", "\\1", Rdata)
    cohort.names = c(cohort.names, cohort_name)
    all.cohort.scores[[cohort_name]] = cohort.data[['score']]
    all.cohort.covms[[cohort_name]] = cohort.data[['covms']]
}
rm(cohort.data)


SMMAT.config = snakemake@params[['SMMAT_config']]
Burden = SMMAT.config[['Burden']]
SKATO = SMMAT.config[['SKATO']]
SKAT = SMMAT.config[['SKAT']]
SMMAT = SMMAT.config[['SMMAT']]
use.minor.allele = SMMAT.config[['use_minor_allele']]
miss.cutoff = SMMAT.config[['miss_cutoff']]
rho = SMMAT.config[['rho']]
method = SMMAT.config[['method']]
MAF.weights.beta = SMMAT.config[['MAF_weights_beta']]
MAF.range = as.numeric(SMMAT.config[['MAF_range']])


subset.groups = fread(snakemake@input[['groups']], header=F)$V1
n.groups = length(subset.groups)

group.info = fread(snakemake@input[['groupfile']], header = FALSE, colClasses = list(
    character = c(1, 2, 4, 5),
    numeric = c(3, 6)
))
setnames(group.info, c("group", "chr", "pos", "ref", "alt", "weight"))

group.info = group.info[group %in% subset.groups]
group.info[,variant.id:=paste(group, chr, format(pos, scientific = FALSE), ref, alt, sep = ":")]
if (any(duplicated(group.info$variant.id))) {
    warning("Duplicated variant ID detected in group.file. Excluding all duplicates, but this should be checked!!")
    group.info <- group.info[!duplicated(variant.id), ]
}

# Below is needed for memory efficient subsetting of the group.info later using variant.ids
setkey(group.info, variant.id)




# Register cores
registerDoMC(cores = snakemake@threads)

# Parallelized loop
results <- foreach(i = 1:n.groups,
                    .combine = rbind,
                    .export = c(
                        "subset.groups",
                        "group.info",
                        "all.cohort.scores",
                        "all.cohort.covms",
                        "MAF.range",
                        "rho",
                        "miss.cutoff",
                        "method",
                        "use.minor.allele",
                        "Burden", "SKATO", "SMMAT", "SKAT")) %dopar% {
    selected.group <- subset.groups[i]

    # Combine summary statistics across cohort score files
    cohort_name = cohort.names[1]
    group.combined.score <- all.cohort.scores[[cohort_name]][group == selected.group, .(variant.id, N, Nmiss, AC, SCORE)]
    
    included.cohorts = NULL
    if (nrow(group.combined.score) > 0) included.cohorts = c(cohort_name)
    for (cohort_name in cohort.names[2:length(cohort.names)]) {
        next.cohort.group.score <- all.cohort.scores[[cohort_name]][group == selected.group, .(variant.id, N, Nmiss, AC, SCORE)]
        if (nrow(next.cohort.group.score) > 0) included.cohorts = c(included.cohorts, cohort_name)
        combined <- merge(group.combined.score, next.cohort.group.score, by = "variant.id", all = TRUE)
        combined[, N := rowSums(.SD, na.rm = TRUE), .SDcols = c("N.x", "N.y")]
        combined[, Nmiss := rowSums(.SD, na.rm = TRUE), .SDcols = c("Nmiss.x", "Nmiss.y")]
        combined[, AC := rowSums(.SD, na.rm = TRUE), .SDcols = c("AC.x", "AC.y")]
        combined[, SCORE := rowSums(.SD, na.rm = TRUE), .SDcols = c("SCORE.x", "SCORE.y")]
        group.combined.score <- combined[, .(variant.id, N, Nmiss, AC, SCORE)]
    }
    included.cohorts.pasted = paste(included.cohorts, collapse=',')
    if (nrow(group.combined.score) == 0) {
        output = data.frame(group = selected.group, n.variants = 0, cohorts = NA)
        if (Burden | SKATO | SMMAT) {
            output$B.score <- NA
            output$B.var <- NA
            output$B.pval <- NA
        }
        if (SKAT | SKATO) output$S.pval <- NA
        if (SKATO) {
            output$O.pval <- NA
            output$O.minp <- NA
            output$O.minp.rho <- NA
        }
        if(SMMAT) output$E.pval <- NA
        return(output)
    }
    
    group.combined.score[, AF := AC / N / 2]
    group.combined.score[, `:=`(
        below.miss.cutoff = (Nmiss / (N + Nmiss) <= miss.cutoff),
        AF.within.range = (MAF.range[1] <= AF & AF <= MAF.range[2]) | (1 - MAF.range[2] <= AF & AF <= 1 - MAF.range[1])
    )]
    group.combined.score <- group.combined.score[below.miss.cutoff & AF.within.range]
    n.variants <- nrow(group.combined.score)
    if (n.variants == 0) {
        # TODO: This is repeated. maybe create a function and just use that.
        output = data.frame(group = selected.group, n.variants = 0, cohorts = NA)
        if (Burden | SKATO | SMMAT) {
            output$B.score <- NA
            output$B.var <- NA
            output$B.pval <- NA
        }
        if (SKAT | SKATO) output$S.pval <- NA
        if (SKATO) {
            output$O.pval <- NA
            output$O.minp <- NA
            output$O.minp.rho <- NA
        }
        if(SMMAT) output$E.pval <- NA
        return(output)
    }
    
    U <- group.combined.score$SCORE
    
    # Combine covariance matrices across cohorts
    V.side <- nrow(group.combined.score)
    V <- matrix(0, V.side, V.side)
    rownames(V) <- group.combined.score$variant.id
    colnames(V) <- group.combined.score$variant.id
    for (cohort_name in cohort.names) {
        cohort.V <- all.cohort.covms[[cohort_name]][[selected.group]]
        IDX <- all.cohort.scores[[cohort_name]][group == selected.group, variant.id]
        if (length(IDX) == 0) next
        colnames(cohort.V) <- IDX
        rownames(cohort.V) <- IDX
        # We first assign the IDX as col/row names for cohort.V, AND THEN we filter out
        # any variants not included in `group.combined.score`.
        # This is because the cohort score file and cohort.V sides follow the same order of variants.    
        IDX2 <- IDX[IDX %in% group.combined.score$variant.id]
        V[IDX2, IDX2] <- V[IDX2, IDX2] + cohort.V[IDX2, IDX2]
    }

    subset.group.info <- group.info[group.combined.score$variant.id] # We can efficiently subset group.info here because of the previous `setkey` call
    group.combined.score.weight = merge(group.combined.score, subset.group.info[, .(variant.id, weight)], by = "variant.id")
    if (use.minor.allele & nrow(group.combined.score.weight[AF>0.05])) {
        group.combined.score.weight[AF>0.05, weight:=-weight]
    }
    group.combined.score.weight[,updated.weight := weight * GMMAT:::MAF.weights.beta.fun(AF, MAF.weights.beta[1], MAF.weights.beta[2])]
    weight = group.combined.score.weight$updated.weight

    U <- U * weight
    V <- t(V*weight)*weight
    output = data.frame(group = selected.group, n.variants = n.variants, cohorts=included.cohorts.pasted)
    if(max(V)-min(V) < sqrt(.Machine$double.eps)) {
        burden.score <- sum(U)
        burden.var <- sum(V)
        burden.pval <- pchisq(burden.score^2/burden.var, df=1, lower.tail=FALSE)
        if (Burden | SKATO | SMMAT) {
            output$B.score <- burden.score
            output$B.var <- burden.var
            output$B.pval <- burden.pval
        }
        if (SKAT | SKATO) output$S.pval <- burden.pval
        if (SKATO) {
            output$O.pval <- burden.pval
            output$O.minp <- burden.pval
            output$O.minp.rho <- 1
        }
        if(SMMAT) output$E.pval <- burden.pval
    } else {
        if(SKATO) {
            re <- GMMAT:::.skato_pval(U = U, V = V, rho = rho, method = method)
            output$B.score <- re$Burden.score
            output$B.var <- re$Burden.var
            output$B.pval <- re$Burden.pval
            output$S.pval <- re$SKAT.pval
            output$O.pval <-re$p
            output$O.minp <- re$minp
            output$O.minp.rho <- re$minp.rho
        } else {
            if(SKAT) output$S.pval <- GMMAT:::.quad_pval(U = U, V = V, method = method)
            if(Burden | SMMAT) {
                output$B.score <- sum(U)
                output$B.var <- sum(V)
                output$B.pval <- pchisq(output$B.score^2/output$B.var, df=1, lower.tail=FALSE)
            }
        }	    
        if(SMMAT) {
            V.rowSums <- rowSums(V)
            U <- U - V.rowSums * output$B.score / output$B.var
            V <- V - tcrossprod(V.rowSums) / output$B.var
            if(mean(abs(V)) < sqrt(.Machine$double.eps)) output$E.pval <- output$B.pval
            else output$E.pval <- tryCatch(pchisq(-2*log(output$B.pval)-2*log(GMMAT:::.quad_pval(U = U, V = V, method = method)), df = 4, lower.tail = FALSE), error = function(e) { output$B.pval })
        }
    }
    print(paste("Completed index:", i, "; group:", selected.group))
    return(output)
}
sink()
# Combine results into final output
results = results[match(subset.groups, results$group), ]

if (any(is.na(results$n.variants))) {
    cat(paste0("Groups which returned NA during the parallel analysis: ", results[is.na(n.variants), group]))
    stop("Some groups failed the analysis!")
}

fwrite(results, file = snakemake@output[["result"]])
