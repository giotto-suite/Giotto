# Bag models ####
#
# A cell-level classifier trained from group-level targets: each group (a
# spot, a bin, a region) holds a bag of cells, the model predicts cell type
# probabilities p_i for every cell, and the losses compare the bag's mean
# prediction pbar with what is known about the group (its deconvolved
# proportions, its expression). Used by runNucleusDeconv() to type H&E nuclei
# from Visium spots; nothing here knows about Visium, images or DWLS.
#
# Three kinds of parts, kept ad hoc for now and shaped like the param objects
# (`method` + `param`) so they can become `analyzeParam` classes later:
#
#   learner  init(d, types) -> theta; forward(theta, X) -> logits;
#            backward(theta, X, dZ) -> gradient for theta (same names)
#   loss     value_grad(P, pbar, ids) -> list(value, d_pbar, d_P), either
#            gradient NULL when the loss does not use it; `weight`; `phase`
#            ("all", or "warm" / "main" for a warm start, see .bag_train())
#   trainer  .bag_train(): Adam over batches of groups, groups split into
#            train / validation / test, the epoch with the lowest validation
#            loss is kept.
#
# The softmax sits between the learner and the losses, so a learner only
# backpropagates through its logits.


# learners ####

#' @keywords internal
#' @noRd
.bag_learner <- function(method = "linear", ...) {
    method <- match.arg(method, choices = "linear")
    switch(method,
        "linear" = .bag_learner_linear(...)
    )
}

# Softmax regression, initialised as torch.nn.Linear (uniform within
# 1 / sqrt(d)).
#' @keywords internal
#' @noRd
.bag_learner_linear <- function() {
    list(
        method = "linear",
        param = list(),
        init = function(d, types) {
            k <- length(types)
            bd <- 1 / sqrt(d)
            list(
                W = matrix(stats::runif(d * k, -bd, bd), d, k,
                    dimnames = list(NULL, types)),
                b = stats::setNames(stats::runif(k, -bd, bd), types)
            )
        },
        forward = function(theta, X) {
            sweep(X %*% theta$W, 2L, theta$b, "+")
        },
        backward = function(theta, X, dZ) {
            list(W = crossprod(X, dZ), b = colSums(dZ))
        }
    )
}


# losses ####

#' @keywords internal
#' @noRd
.bag_loss <- function(method = c("proportion", "dwls", "confidence"),
    ..., weight = 1, phase = c("all", "warm", "main")) {
    method <- match.arg(method)
    phase <- match.arg(phase)
    l <- switch(method,
        "proportion" = .bag_loss_proportion(...),
        "dwls" = .bag_loss_dwls(...),
        "confidence" = .bag_loss_confidence(...)
    )
    l$weight <- weight
    l$phase <- phase
    l
}

# Mean squared difference between pbar and given proportions (groups x
# types, rows named by group). HEDeST's proportion loss.
.bag_loss_proportion <- function(target) {
    target <- as.matrix(target)
    if (is.null(rownames(target)) || is.null(colnames(target))) {
        stop("[bag model] proportion targets need group IDs as row names ",
            "and cell types as column names.", call. = FALSE)
    }
    list(
        method = "proportion",
        param = list(target = target),
        groups = rownames(target),
        value_grad = function(P, pbar, ids) {
            d <- pbar - target[ids, colnames(pbar), drop = FALSE]
            list(value = mean(d^2), d_pbar = 2 * d / length(d), d_P = NULL)
        }
    )
}

# spatialDWLS's dampened WLS of each group's expression on pbar, from
# spatialDWLS::dwls_design().
.bag_loss_dwls <- function(design) {
    .dwls_engine_check()
    if (!exists("dwls_expression_loss", envir = asNamespace("spatialDWLS"),
            inherits = FALSE)) {
        stop("[bag model] the DWLS expression loss needs spatialDWLS >= ",
            "0.1.2: remotes::install_github(\"drieslab/spatialDWLS\")",
            call. = FALSE)
    }
    list(
        method = "dwls",
        param = list(design = design),
        groups = colnames(design$y),
        value_grad = function(P, pbar, ids) {
            l <- spatialDWLS::dwls_expression_loss(design, t(pbar))
            list(value = l$value, d_pbar = t(l$grad), d_P = NULL)
        }
    )
}

# -mean(log max p_i): pushes each cell towards one type. HEDeST's tie-breaker.
.bag_loss_confidence <- function() {
    list(
        method = "confidence",
        param = list(),
        groups = NULL,
        value_grad = function(P, pbar, ids) {
            am <- max.col(P, ties.method = "first")
            idx <- cbind(seq_len(nrow(P)), am)
            pm <- P[idx]
            d_P <- matrix(0, nrow(P), ncol(P))
            d_P[idx] <- -1 / (nrow(P) * pm)
            list(value = -mean(log(pm)), d_pbar = NULL, d_P = d_P)
        }
    )
}


# training ####

.bag_softmax <- function(Z) {
    Z <- Z - apply(Z, 1L, max)
    E <- exp(Z)
    E / rowSums(E)
}

# Loss over the groups `ids` and, with grad = TRUE, its gradient for theta.
# `cells` maps every group to row indices of X.
.bag_batch <- function(theta, learner, X, cells, ids, losses, grad = TRUE) {
    idx <- unlist(cells[ids], use.names = FALSE)
    ns <- lengths(cells[ids])
    bag <- rep.int(seq_along(ids), ns)
    Xb <- X[idx, , drop = FALSE]
    P <- .bag_softmax(learner$forward(theta, Xb))
    pbar <- rowsum(P, bag, reorder = FALSE) / ns
    dimnames(pbar) <- list(ids, colnames(P))
    value <- 0
    d_pbar <- matrix(0, nrow(pbar), ncol(pbar))
    d_P <- matrix(0, nrow(P), ncol(P))
    for (l in losses) {
        r <- l$value_grad(P, pbar, ids)
        value <- value + l$weight * r$value
        if (!is.null(r$d_pbar)) d_pbar <- d_pbar + l$weight * r$d_pbar
        if (!is.null(r$d_P)) d_P <- d_P + l$weight * r$d_P
    }
    if (!isTRUE(grad)) return(value)
    d_P <- d_P + d_pbar[bag, , drop = FALSE] / ns[bag]
    dZ <- P * (d_P - rowSums(d_P * P))
    list(value = value, grad = learner$backward(theta, Xb, dZ))
}

# Train `learner` on X (cells x features, row names cell IDs) from `groups`
# (data.frame of group + cell_ID, or a named list of cell IDs per group).
# Only groups that every loss has a target for, and that hold at least one
# cell of X, are used. `warm` > 0: the first `warm` epochs use the losses
# with phase "warm" (and "all"), the rest those with phase "main" (and
# "all"); the kept epoch is chosen within the main phase.
#' @keywords internal
#' @noRd
.bag_train <- function(X, groups, losses, learner = .bag_learner("linear"),
    types, epochs = 100L, lr = 1e-3, batch_size = 64L, wd = 0,
    split = c(train = 0.8, val = 0.1), warm = 0L, seed = NULL,
    verbose = NULL) {
    if (is.null(rownames(X))) {
        stop("[bag model] X needs cell IDs as row names.", call. = FALSE)
    }
    if (!is.null(seed)) GiottoUtils::local_seed(seed)
    cells <- .bag_cells(groups, rownames(X))
    for (l in losses) {
        if (!is.null(l$groups)) cells <- cells[names(cells) %in% l$groups]
    }
    if (length(cells) < 10L) {
        stop("[bag model] fewer than 10 groups hold cells and have targets ",
            "for every loss.", call. = FALSE)
    }
    ids <- sample(names(cells))
    n_tr <- floor(split[[1]] * length(ids))
    n_va <- max(1L, floor(split[[2]] * length(ids)))
    tr <- ids[seq_len(n_tr)]
    va <- ids[n_tr + seq_len(n_va)]
    active <- function(phase) {
        Filter(function(l) l$phase %in% c("all", phase), losses)
    }

    theta <- learner$init(ncol(X), types)
    mom <- lapply(theta, function(x) x * 0)
    vel <- mom
    step <- 0L
    best <- Inf
    best_theta <- theta
    best_epoch <- NA_integer_
    curve <- numeric(epochs)
    for (ep in seq_len(epochs)) {
        phase <- if (warm > 0L && ep <= warm) "warm" else "main"
        if (warm > 0L && ep == warm + 1L) best <- Inf
        L <- active(phase)
        o <- sample(tr)
        for (s0 in seq(1L, length(o), by = batch_size)) {
            b <- o[s0:min(s0 + batch_size - 1L, length(o))]
            g <- .bag_batch(theta, learner, X, cells, b, L)$grad
            step <- step + 1L
            for (p in names(theta)) {
                gp <- g[[p]] + wd * theta[[p]]
                mom[[p]] <- 0.9 * mom[[p]] + 0.1 * gp
                vel[[p]] <- 0.999 * vel[[p]] + 0.001 * gp^2
                theta[[p]] <- theta[[p]] - lr *
                    (mom[[p]] / (1 - 0.9^step)) /
                    (sqrt(vel[[p]] / (1 - 0.999^step)) + 1e-8)
            }
        }
        curve[ep] <- .bag_batch(theta, learner, X, cells, va, L, grad = FALSE)
        if (curve[ep] < best) {
            best <- curve[ep]
            best_theta <- theta
            best_epoch <- ep
        }
    }
    vmsg(.v = verbose, sprintf(
        "[bag model] %d groups (%d train, %d validation), best epoch %d",
        length(ids), length(tr), length(va), best_epoch))
    list(
        learner = learner, theta = best_theta, types = types,
        best_epoch = best_epoch, val_loss = curve,
        groups = list(train = tr, val = va,
            test = setdiff(ids, c(tr, va)))
    )
}

# Cell type probabilities for the rows of X (cells x features), in chunks.
#' @keywords internal
#' @noRd
.bag_predict <- function(model, X, chunk_size = 50000L) {
    starts <- seq(1L, max(nrow(X), 1L), by = chunk_size)
    P <- do.call(rbind, lapply(starts, function(s) {
        i <- s:min(s + chunk_size - 1L, nrow(X))
        .bag_softmax(model$learner$forward(model$theta,
            X[i, , drop = FALSE]))
    }))
    dimnames(P) <- list(rownames(X), model$types)
    P
}

# Groups as a named list of row indices into the cells `ids`; cells not in
# `ids` are dropped, and so are groups left empty.
.bag_cells <- function(groups, ids) {
    if (is.data.frame(groups)) {
        if (!all(c("group", "cell_ID") %in% colnames(groups))) {
            stop("[bag model] a groups table needs columns `group` and ",
                "`cell_ID`.", call. = FALSE)
        }
        groups <- split(as.character(groups$cell_ID),
            as.character(groups$group))
    }
    out <- lapply(groups, function(g) {
        i <- match(g, ids)
        i[!is.na(i)]
    })
    out[lengths(out) > 0L]
}
