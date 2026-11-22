# Regenerate numerical fixtures with an independently installed ade4.
# Run from the package root: Rscript tools/make_reference.R
library(ade4)
options(digits = 17)
out <- file.path('test', 'fixtures')
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_matrix <- function(x, name) {
  write.table(as.matrix(x), file.path(out, paste0(name, '.tsv')),
              sep = '\t', row.names = FALSE, col.names = FALSE, quote = FALSE, na = 'NaN')
}
interactive_call <- function(fun, answers, ...) {
  env <- new.env(parent = environment(fun))
  env$readLines <- function(...) {
    answer <- as.character(answers[1])
    answers <<- answers[-1]
    answer
  }
  environment(fun) <- env
  invisible(capture.output(result <- fun(..., scann = TRUE)))
  result
}
export_case <- function(name, tables, types, answers = NULL, option = 'scaledBYrange') {
  tab <- ktab.list.df(tables)
  call <- function(fun, ...) {
    if (is.null(answers)) fun(tab, types, option = option, ...)
    else interactive_call(fun, answers, tab, types, option = option, ...)
  }
  write_matrix(call(dist.ktab), paste0(name, '_global'))
  per_trait <- call(ldist.ktab)
  write_matrix(do.call(cbind, lapply(per_trait, as.vector)), paste0(name, '_traits'))
  for (squared in c(TRUE, FALSE)) {
    correlations <- suppressWarnings(call(kdist.cor, squared = squared))
    for (field in names(correlations))
      write_matrix(correlations[[field]], paste0(name, '_', if(squared) 'sq_' else 'raw_', field))
  }
}

Q <- data.frame(q1 = c(1, 2, 4, 8, 5, 3), q2 = c(0, 1, 7, 2, 4, 6))
O <- data.frame(o1 = c(1, 1, 2, 3, 3, 4), o2 = c(0, 2, 2, 2, 4, 4))
N <- data.frame(n1 = c('a','b','a','c','b','a'), n2 = c('x','x','y','y','z','z'))
D <- data.frame(d1 = c(0,1,0,1,0,1), d2 = c(1,1,0,0,1,0))
F <- data.frame(f1 = c(1,2,0,1,4,2), f2 = c(0,1,3,1,1,0), f3 = c(2,0,1,1,0,1),
                f4 = c(1,3,0,2,1,1), f5 = c(1,0,1,1,4,2))
B <- data.frame(b1 = c(1,0,1,1,0,1), b2 = c(0,1,1,0,1,0), b3 = c(1,1,0,1,1,1),
                b4 = c(1,0,1,1,0,1), b5 = c(0,1,1,1,1,0))
C <- data.frame(c1 = c(1,12,6,7,3,10), c2 = c(1,5,3,4,2,1))
raw <- list(Q=Q, O=O, N=N, D=D, F=F, B=B, C=C)
for (name in names(raw)) write_matrix(raw[[name]], paste0('input_', name))
prepare <- function(x) {
  x$F <- prep.fuzzy(x$F, c(3,2), labels = c('F1','F2'))
  x$B <- prep.binary(x$B, c(3,2), labels = c('B1','B2'))
  x$C <- prep.circular(x$C, c(1,1), c(12,5))
  x
}
export_case('mixed', prepare(raw), names(raw))
missing <- raw
missing$Q[2,1] <- NA
missing$O[5,2] <- NA
missing$N[3,1] <- NA
missing$D[1,2] <- NA
missing$F[4,1:3] <- NA
missing$B[6,4:5] <- NA
missing$C[2,1] <- NA
for (name in names(missing)) write_matrix(missing[[name]], paste0('missing_', name))
export_case('missing', prepare(missing), names(raw))
export_case('noscale', prepare(raw), names(raw), option = 'noscale')
export_case('q_sd', list(Q), 'Q', option = 'scaledBYsd')
# Interactive answers: F, B, O, Q (when that type is present).
export_case('manhattan', prepare(raw), names(raw), c(2,1,1,2))
export_case('podani', prepare(raw), names(raw), c(2,1,2,1))
for (f in c(1,3,4))
  export_case(paste0('fuzzy',f), list(prep.fuzzy(F,c(3,2))), 'F', f)
Bstable <- as.data.frame(rbind(c(1,0,0,0), c(0,1,1,0), c(1,0,1,0),
                              c(0,1,0,1), c(1,1,0,1), c(0,0,0,1)))
write_matrix(Bstable, 'binary_input')
for (b in 1:10)
  export_case(paste0('binary',b), list(prep.binary(Bstable,4)), 'B', b)

data(woangers)
t <- woangers$traits
write_matrix(t, 'woangers_input')
wo <- list(t[,c(1,2,7,8)], prep.circular(t[3],1,12), t[,4:6],
           prep.fuzzy(t[,9:19],c(3,3,5),labels=c('mp','pe','di')), t[,20:21])
export_case('woangers', wo, c('N','C','O','F','Q'))
writeLines(paste0('ade4 ', packageVersion('ade4'), '; R ', getRversion()),
           file.path(out, 'VERSION'))
cat('Exported fixtures from ade4', as.character(packageVersion('ade4')), '\n')
