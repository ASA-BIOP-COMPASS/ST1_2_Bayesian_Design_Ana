#' Convert a prespecified fraction of total control information into borrowing ESS
#' @param n_control Number of randomized concurrent control patients.
#' @param borrowing_prop Fraction of total nominal control information contributed by historical data, in [0,1).
#' @return Numeric historical ESS, which may be fractional.
#' @export
borrowing_ess <- function(n_control, borrowing_prop = 0) {
  .nonnegative(n_control, 'n_control')
  if (n_control != floor(n_control)) stop('n_control must be an integer.', call. = FALSE)
  if (!is.numeric(borrowing_prop) || length(borrowing_prop) != 1L ||
      !is.finite(borrowing_prop) || borrowing_prop < 0 || borrowing_prop >= 1)
    stop('borrowing_prop must be within [0,1).', call. = FALSE)
  if (n_control == 0 && borrowing_prop > 0)
    stop('A positive borrowing proportion needs concurrent controls. For historical-only control use a fixed prior in analyze_binary().', call. = FALSE)
  n_control * borrowing_prop / (1 - borrowing_prop)
}

#' Resolve exactly one borrowing specification
#' @keywords internal
.resolve_borrow <- function(n_control, borrowed_n, borrowing_prop) {
  if (!is.null(borrowed_n) && !is.null(borrowing_prop))
    stop('Specify borrowed_n OR borrowing_prop, not both.', call. = FALSE)
  if (!is.null(borrowed_n)) {
    .nonnegative(borrowed_n, 'borrowed_n')
    return(borrowed_n)
  }
  borrowing_ess(n_control, if (is.null(borrowing_prop)) 0 else borrowing_prop)
}

#' Calculate overall power for one or more treatment arms
#' @param n Named recruited sample sizes, placebo first.
#' @param response_rates Named anticipated response rates.
#' @param borrowed_n Fixed historical effective sample size, OR use borrowing_prop.
#' @param borrowing_prop Borrowed ESS/(borrowed ESS + recruited placebo n).
#' @param historical_rate Assumed historical placebo rate, default concurrent rate.
#' @param overall One or all doses must succeed.
#' @param multiplicity Multiplicity method; default 'hochberg'.
#' @param seed Reproducible simulation seed.
#' @param base_control_prior Baseline Beta prior shapes for placebo;
#'   default c(0.5, 0.5), before any explicit historical borrowing.
#' @param treatment_prior Common Beta prior shapes for each active dose;
#'   default c(0.5, 0.5).
#' @param ... Passed to simulate_power.
#' @export
design_power <- function(n, response_rates, borrowed_n = NULL,
                         borrowing_prop = NULL, historical_rate = NULL,
                         overall = c('any','all'), multiplicity = 'hochberg',
                         seed = 123, base_control_prior = c(0.5, 0.5),
                         treatment_prior = c(0.5, 0.5), ...) {
  overall <- match.arg(overall)
  if (is.null(names(n)) || !identical(names(n), names(response_rates)))
    stop('n and response_rates must have identical ordered names.', call. = FALSE)
  ess <- .resolve_borrow(n[1], borrowed_n, borrowing_prop)
  fit <- simulate_power(n, response_rates, borrow_ess = ess,
                        historical_rate = historical_rate, multiplicity = multiplicity,
                        seed = seed, base_control_prior = base_control_prior,
                        treatment_prior = treatment_prior, ...)
  fit$overall_power <- if (overall == 'any') fit$any_power else fit$all_power
  fit$overall_definition <- overall
  fit$borrowed_ess <- ess
  fit$borrowing_prop <- if (n[1] + ess > 0) ess / (n[1] + ess) else 0
  fit
}

#' Power-based sample-size search using exact integer allocation blocks
#'
#' Candidate totals not divisible by allocation_unit(ratios) are skipped,
#' rather than rounding individual arm sizes. To generate only valid totals,
#' use candidate_totals(ratios, min_total, max_total).
#'
#' With search='two_stage', evaluate every initial_step-th feasible candidate
#' using n_sim_initial trials. Find the earliest preliminary candidate reaching
#' target_power (or the closest candidate if none reaches it), and refine a
#' window of refine_span feasible candidates either side using n_sim trials.
#' Final sample-size selection uses ONLY final-stage simulation results.
#' `results` combines early and refined estimates, replacing overlapping
#' preliminary estimates with final estimates; `initial_results` and
#' `final_results` allow auditing the two separate stages.
#'
#' @param total_n Candidate randomized totals. Non-feasible totals are omitted.
#' @param ratios Named positive integer allocation weights, placebo first.
#' @param response_rates Named anticipated response probabilities, same names/order.
#' @param target_power Desired overall power, by default 0.8.
#' @param borrowed_n Fixed borrowed historical ESS; mutually exclusive with borrowing_prop.
#' @param borrowing_prop Borrowed ESS / (borrowed ESS + concurrent placebo n).
#' @param historical_ess Optional conservative MAP ESS reference, never used for design optimization.
#' @param plot Draw the stored power results as a curve plus transposed table.
#' @param overall 'any' (at least one dose), or 'all' (every dose).
#' @param multiplicity Posterior-tail score adjustment: none, bonferroni, holm, hochberg (default), BH.
#' @param seed Simulation seed; default 123.
#' @param base_control_prior Baseline placebo Beta shapes before borrowing, default c(0.5,0.5).
#' @param treatment_prior Common treatment Beta shapes, default c(0.5,0.5).
#' @param search 'grid' (all valid candidates at full n_sim) or 'two_stage' (coarse/final).
#' @param n_sim Trial simulations per grid candidate or FINAL two-stage candidate.
#' @param n_sim_initial Number of trial simulations per coarse-stage candidate.
#' @param initial_step Every initial_step-th feasible candidate for coarse stage.
#' @param refine_span Number of feasible candidates either side of coarse anchor.
#'   By default equals initial_step; refinements include an earlier coarse point.
#' @param table_max_cols Maximum columns displayed in the embedded plot table.
#' @param ... Additional simulate_power arguments (n_post, margin, alpha, historical_rate).
#' @return A bda_sample_size object with results, initial_results, final_results,
#'   minimum_candidate, and note. Monte Carlo SE is estimated for each row.
#' @examples
#' ratio <- c(Placebo=1, Low=1, Mid=1, High=1)
#' grid <- candidate_totals(ratio, 160, 240)
#' # ss <- design_sample_size(grid, ratio, c(Placebo=.3, Low=.4, Mid=.5, High=.6),
#' #   search='two_stage', n_sim_initial=200, n_sim=2000, plot=TRUE)
#' @export
design_sample_size <- function(total_n, ratios, response_rates,
                               target_power=.8, borrowed_n=NULL,
                               borrowing_prop=NULL, historical_ess=NULL,
                               plot=FALSE, overall=c('any','all'),
                               multiplicity='hochberg', seed=123,
                               base_control_prior=c(.5,.5),
                               treatment_prior=c(.5,.5),
                               search=c('grid','two_stage'), n_sim=1000,
                               n_sim_initial=500, initial_step=5L,
                               refine_span=NULL, table_max_cols=12L, ...) {
  overall <- match.arg(overall)
  search <- match.arg(search)
  .prob(target_power,'target_power')
  if (!identical(names(ratios),names(response_rates)))
    stop('ratios and response_rates must have identical arm names in the same order.',call.=FALSE)
  block <- allocation_unit(ratios)
  .prob(response_rates,'response_rates')
  if (!is.numeric(total_n) || !length(total_n) || any(!is.finite(total_n)) ||
      any(total_n < 0) || any(total_n != floor(total_n)))
    stop('total_n must contain nonnegative integers.',call.=FALSE)
  if (!is.null(historical_ess)) .nonnegative(historical_ess,'historical_ess')
  # Reject two competing borrowing methods before doing any simulations.
  .resolve_borrow(block, borrowed_n, borrowing_prop)
  .validate_count <- function(x, label) {
    if (!is.numeric(x) || length(x)!=1L || !is.finite(x) || x < 1 || x != floor(x))
      stop(paste0(label,' must be a positive integer.'),call.=FALSE)
  }
  for (arg in c('n_sim','n_sim_initial','initial_step','table_max_cols'))
    .validate_count(get(arg),arg)
  if (is.null(refine_span)) refine_span <- initial_step
  .validate_count(refine_span,'refine_span')
  feasible <- sort(unique(as.integer(total_n[total_n >= block & total_n %% block == 0])))
  omitted <- sort(unique(total_n[!total_n %in% feasible]))
  if (length(omitted))
    warning(sprintf('Skipped %d total sample sizes incompatible with the exact ratio (block=%d). Use candidate_totals().',
                    length(omitted), block), call.=FALSE)
  if (!length(feasible))
    stop(sprintf('No feasible totals. This ratio requires multiples of %d. Use candidate_totals().', block),call.=FALSE)
  dots <- list(...)
  if (any(names(dots) %in% c('n','response_rates','borrowed_n','borrowing_prop',
                            'overall','multiplicity','seed','base_control_prior',
                            'treatment_prior','n_sim')))
    stop('Do not repeat named design arguments in ...',call.=FALSE)
  evaluate <- function(candidate, simulations, stage) {
    rows <- lapply(candidate, function(total) {
      nn <- allocate_arms(total,ratios)
      fit <- do.call(design_power,c(list(n=nn, response_rates=response_rates,
        borrowed_n=borrowed_n, borrowing_prop=borrowing_prop,
        overall=overall, multiplicity=multiplicity, seed=seed,
        base_control_prior=base_control_prior, treatment_prior=treatment_prior,
        n_sim=simulations),dots))
      data.frame(total_n=as.integer(total), overall_power=fit$overall_power,
        mc_se=sqrt(fit$overall_power*(1-fit$overall_power)/simulations),
        borrowed_ess=fit$borrowed_ess, borrowing_prop=fit$borrowing_prop,
        as.list(nn), n_sim=as.integer(simulations), stage=stage,
        check.names=FALSE)
    })
    tab <- do.call(rbind,rows)
    rownames(tab) <- NULL
    tab$meets_target <- tab$overall_power >= target_power
    if (!is.null(historical_ess)) {
      tab$historical_ess_reference <- historical_ess
      tab$exceeds_historical_ess <- tab$borrowed_ess > historical_ess + 1e-8
    }
    tab
  }
  if (search=='grid') {
    initial <- NULL
    refined <- evaluate(feasible,n_sim,'full')
    tab <- refined
    anchor_n <- NA_integer_
  } else {
    coarse_indices <- unique(c(seq.int(1L,length(feasible),by=initial_step), length(feasible)))
    initial <- evaluate(feasible[coarse_indices],n_sim_initial,'initial')
    passing <- which(initial$meets_target)
    if (length(passing)) {
      anchor_n <- initial$total_n[passing[1L]]
    } else {
      # Estimate most promising area while also covering the upper grid edge.
      anchor_n <- initial$total_n[which.max(initial$overall_power)]
    }
    anchor_pos <- match(anchor_n,feasible)
    region <- seq.int(max(1L,anchor_pos-refine_span),
                      min(length(feasible),anchor_pos+refine_span))
    refined <- evaluate(feasible[region],n_sim,'final')
    tab <- rbind(initial[!initial$total_n %in% refined$total_n,,drop=FALSE],refined)
    tab <- tab[order(tab$total_n),,drop=FALSE]
    rownames(tab) <- NULL
  }
  # Only the high precision trials may drive final selection.
  passing_final <- which(refined$meets_target)
  chosen <- if (length(passing_final)) refined[passing_final[1L],,drop=FALSE] else NULL
  note <- NULL
  if (is.null(chosen)) {
    note <- 'No FINAL-stage candidate reached the target; expand the grid or increase the refinement window.'
  } else if (!is.null(historical_ess) && chosen$exceeds_historical_ess) {
    note <- sprintf('NOTE: Borrowed ESS (%.2f) exceeds historical reference ESS (%.2f). This is a warning only, not a design constraint.',
                    chosen$borrowed_ess,historical_ess)
  }
  out <- list(results=tab, initial_results=initial, final_results=refined,
    minimum_candidate=chosen, target_power=target_power,
    overall_definition=overall, multiplicity=multiplicity, seed=seed,
    ratios=ratios, allocation_block=block, search=search,
    initial_step=initial_step, refine_span=refine_span, anchor_n=anchor_n,
    n_sim_final=n_sim, n_sim_initial=if(search=='two_stage')n_sim_initial else NULL,
    table_max_cols=table_max_cols, historical_ess_reference=historical_ess,
    base_control_prior=base_control_prior, treatment_prior=treatment_prior,
    note=note)
  class(out) <- 'bda_sample_size'
  if (plot) graphics::plot(out)
  out
}

# Pick a readable set of columns near the target rather than rendering
# 60+ microscopic table columns. The underlying results always remain complete.
.focus_table_columns <- function(totals, selected, maximum) {
  if (length(totals) <= maximum) return(seq_along(totals))
  key <- if (length(selected) && any(is.finite(selected))) selected[1L] else
    totals[ceiling(length(totals)/2)]
  anchor <- which.min(abs(totals-key))
  lo <- max(1L, min(anchor - floor(maximum / 2), length(totals)-maximum+1L))
  seq.int(lo, lo+maximum-1L)
}

# Render one transposed table with rows=parameters and columns=total N.
.plot_bda_table <- function(headers, row_labels, cells, highlights=NULL,
                            row_colors=NULL) {
  nc <- length(headers); nr <- length(row_labels)
  if (is.null(highlights)) highlights <- matrix(FALSE,nrow=nr,ncol=nc)
  if (is.null(row_colors)) row_colors <- rep('black',nr)
  if (!identical(dim(cells),c(nr,nc)))
    stop('Plot table dimensions are inconsistent.',call.=FALSE)
  label_width <- max(3.5, max(nchar(row_labels)) * .22)
  graphics::plot.new()
  graphics::plot.window(xlim=c(-label_width,nc+.6),ylim=c(0,nr+1),xaxs='i',yaxs='i')
  font <- max(.45,min(.85, 10 / max(10,nc)))
  for (j in seq_len(nc)) {
    graphics::rect(j-.5,nr,j+.5,nr+1,col='#E8ECF1',border='#D3D8DD')
    graphics::text(j,nr+.5,headers[j],cex=font,font=2)
    for(i in seq_len(nr)) {
      yy <- nr-i
      bg <- if(highlights[i,j]) '#DCEFD9' else
        if(i %% 2L == 0) '#F7F8F9' else 'white'
      graphics::rect(j-.5,yy,j+.5,yy+1,col=bg,border='#D3D8DD')
      graphics::text(j,yy+.5,cells[i,j],cex=font,
                     font=if(highlights[i,j]) 2 else 1)
    }
  }
  for (i in seq_len(nr)) graphics::text(.25,nr-i+.5,row_labels[i],
    adj=1,cex=font,col=row_colors[i])
  invisible(NULL)
}

#' Draw a power curve and a transposed allocation table from existing results
#' @param x Object returned by design_sample_size().
#' @param ... Unused.
#' @export
plot.bda_sample_size <- function(x, ...) {
  tab <- x$results
  arms <- names(x$ratios)
  selected_n <- if(is.null(x$minimum_candidate)) numeric() else x$minimum_candidate$total_n
  op <- graphics::par(no.readonly=TRUE)
  on.exit({graphics::layout(1);graphics::par(op)},add=TRUE)
  graphics::layout(matrix(c(1,2),2,1),heights=c(3.2,max(2,.35*(length(arms)+4))))
  graphics::par(mar=c(4.8,4.2,3,1.0))
  graphics::plot(tab$total_n,tab$overall_power,type='b',pch=16,ylim=c(0,1),
    xlab='Total randomized sample size',ylab='Overall power',
    main=paste('Overall power:',x$overall_definition,'dose success'))
  graphics::abline(h=x$target_power,lty=2,col='gray40')
  if (x$search == 'two_stage') {
    preliminary <- tab$stage == 'initial'
    if (any(preliminary)) graphics::points(tab$total_n[preliminary],
      tab$overall_power[preliminary],pch=1,col='gray45',cex=1.05)
    graphics::legend('bottomright',legend=c('Final / high precision','Initial / screening'),
      pch=c(16,1),col=c('black','gray45'),bty='n',cex=.70)
  }
  if (length(selected_n)) {
    graphics::abline(v=selected_n,lty=3,col='#17633B')
    graphics::axis(1,at=selected_n,labels=paste0('N=',selected_n),
      line=1.1,col.axis='#17633B',font=2,cex.axis=.9)
    graphics::points(selected_n,x$minimum_candidate$overall_power,
      pch=21,bg='#DCEFD9',col='#17633B',cex=1.75,lwd=2)
  }
  cols <- .focus_table_columns(tab$total_n,selected_n,x$table_max_cols)
  z <- tab[cols,,drop=FALSE]
  labels <- c('Total N',arms,'Borrowed ESS','Overall power')
  cell <- rbind(as.character(z$total_n),
    do.call(rbind,lapply(arms,function(arm) as.character(z[[arm]]))),
    sprintf('%.1f',z$borrowed_ess),sprintf('%.1f%%',100*z$overall_power))
  hl <- matrix(FALSE,length(labels),nrow(z))
  if(length(selected_n)) hl[,which(z$total_n==selected_n)] <- TRUE
  graphics::par(mar=c(1.2,1.1,1,1.0))
  .plot_bda_table(as.character(z$total_n),labels,cell,hl)
  if(!is.null(x$note)) graphics::mtext(x$note,side=1,line=.12,cex=.51,col='#8A3C22')
  invisible(x)
}

#' Compare power across joint dose-specific risk-difference scenarios
#' @param delta_matrix Matrix, one row per scenario and one column per active dose.
#' @param total_n Candidate totals; only exact integer-ratio allocations are evaluated.
#' @param ratios Named positive integer allocation weights, placebo first.
#' @param placebo_rate Anticipated concurrent placebo response probability.
#' @param borrowed_n Fixed borrowed ESS, alternative to borrowing_prop.
#' @param borrowing_prop Borrowed / (borrowed + concurrent placebo n).
#' @param historical_ess Optional historical MAP ESS reference for warnings only.
#' @param plot If TRUE, render curves and matching transposed allocation table.
#' @param overall 'any' or 'all' dose success, default 'any'.
#' @param multiplicity Posterior tail-score method, default 'hochberg'.
#' @param seed Simulation seed.
#' @param base_control_prior Beta prior shapes before borrowing, default c(.5,.5).
#' @param treatment_prior Beta shapes for active doses, default c(.5,.5).
#' @param table_max_cols Number of displayed columns in transposed allocation table.
#' @param ... Additional design_sample_size() options, including search,
#'   n_sim, n_sim_initial, initial_step, refine_span, n_post, alpha, and margin.
#' @return List with results, scenarios and legend labels. No extra simulations for drawing.
#' @export
plot_power <- function(delta_matrix,total_n,ratios,placebo_rate,
                       borrowed_n=NULL,borrowing_prop=NULL,
                       historical_ess=NULL,plot=TRUE,
                       overall=c('any','all'),multiplicity='hochberg',seed=123,
                       base_control_prior=c(.5,.5),treatment_prior=c(.5,.5),
                       table_max_cols=12L,...) {
  overall <- match.arg(overall)
  .prob(placebo_rate,'placebo_rate')
  allocation_unit(ratios)
  if (is.vector(delta_matrix) && is.numeric(delta_matrix))
    delta_matrix <- matrix(delta_matrix,nrow=1)
  if (!is.matrix(delta_matrix) || !is.numeric(delta_matrix) ||
      ncol(delta_matrix)!=length(ratios)-1L || any(!is.finite(delta_matrix)))
    stop('delta_matrix must have one column per active treatment.',call.=FALSE)
  if (!is.null(colnames(delta_matrix)) &&
      !identical(colnames(delta_matrix),names(ratios)[-1]))
    stop('delta_matrix column names must match treatments in order.',call.=FALSE)
  if (!is.numeric(table_max_cols) || length(table_max_cols)!=1L ||
      !is.finite(table_max_cols) || table_max_cols < 1 ||
      table_max_cols != floor(table_max_cols))
    stop('table_max_cols must be a positive integer.',call.=FALSE)
  curves <- vector('list',nrow(delta_matrix))
  studies <- vector('list',nrow(delta_matrix))
  labels <- character(nrow(delta_matrix))
  for (i in seq_len(nrow(delta_matrix))) {
    rr <- setNames(c(placebo_rate,placebo_rate+delta_matrix[i,]),names(ratios))
    .prob(rr,'derived response rates')
    studies[[i]] <- design_sample_size(total_n,ratios,rr,
      borrowed_n=borrowed_n,borrowing_prop=borrowing_prop,
      historical_ess=historical_ess,plot=FALSE,overall=overall,
      multiplicity=multiplicity,seed=seed,
      base_control_prior=base_control_prior,treatment_prior=treatment_prior,
      table_max_cols=table_max_cols,...)
    labels[i] <- paste(sprintf('%s=%+.2f',names(ratios)[-1],delta_matrix[i,]),collapse=', ')
    curves[[i]] <- transform(studies[[i]]$results,scenario=labels[i])
  }
  tab <- do.call(rbind,curves)
  rownames(tab) <- NULL
  if (plot) {
    op <- graphics::par(no.readonly=TRUE)
    on.exit({graphics::layout(1);graphics::par(op)},add=TRUE)
    cols <- grDevices::hcl.colors(length(curves),'Dark 3')
    totals <- sort(unique(tab$total_n))
    selections <- vapply(studies,function(z)
      if(is.null(z$minimum_candidate)) NA_real_ else z$minimum_candidate$total_n,
      numeric(1))
    graphics::layout(matrix(c(1,2),2,1),
      heights=c(3.3,max(2,.35*(length(ratios)+length(curves)+3))))
    graphics::par(mar=c(4.8,4.3,3,1))
    graphics::plot(range(totals),c(0,1),type='n',
      xlab='Total randomized sample size',ylab='Overall power',
      main=paste('Overall power across scenarios:',overall,'dose success'))
    graphics::abline(h=studies[[1]]$target_power,lty=2,col='gray40')
    for (i in seq_along(curves)) {
      graphics::lines(curves[[i]]$total_n,curves[[i]]$overall_power,
        type='b',pch=16,col=cols[i],lwd=2)
      if(is.finite(selections[i])) {
        z <- studies[[i]]$minimum_candidate
        graphics::points(z$total_n,z$overall_power,pch=21,bg='white',
                         col=cols[i],cex=1.7,lwd=2)
      }
    }
    mark_x <- sort(unique(selections[is.finite(selections)]))
    if(length(mark_x)) graphics::axis(1,at=mark_x,
      labels=paste0('N=',mark_x),line=1.1,col.axis='#17633B',font=2,cex.axis=.76)
    graphics::legend('bottomright',legend=labels,col=cols,lty=1,pch=16,
      bty='n',cex=max(.45,min(.75,6/length(labels))))
    # Each scenario may have a different two-stage refined grid. Join by
    # total N, leaving '—' where one scenario has not been simulated.
    shown <- .focus_table_columns(totals,mark_x,table_max_cols)
    displayed <- totals[shown]
    arm_mat <- t(vapply(displayed,function(total)
      allocate_arms(total,ratios),integer(length(ratios))))
    colnames(arm_mat) <- names(ratios)
    borrowed <- vapply(displayed,function(total) {
      nn <- allocate_arms(total,ratios)
      .resolve_borrow(nn[1],borrowed_n,borrowing_prop)
    },numeric(1))
    powers <- do.call(rbind,lapply(curves,function(curve) {
      vals <- curve$overall_power[match(displayed,curve$total_n)]
      ifelse(is.na(vals),'—',sprintf('%.1f%%',100*vals))
    }))
    row_labs <- c('Total N',names(ratios),'Borrowed ESS',paste0('Power: ',labels))
    cell <- rbind(as.character(displayed),t(arm_mat),
      sprintf('%.1f',borrowed),powers)
    hl <- matrix(FALSE,nrow=length(row_labs),ncol=length(displayed))
    for (i in seq_along(selections)) if(is.finite(selections[i])) {
      hit <- which(displayed==selections[i])
      if(length(hit)) hl[length(ratios)+2L+i,hit] <- TRUE
    }
    graphics::par(mar=c(1.2,1.1,1,1))
    .plot_bda_table(as.character(displayed),row_labs,cell,hl,
      row_colors=c(rep('black',length(ratios)+2L),cols))
  }
  ans <- list(results=tab,scenarios=studies,legend=labels)
  if (plot) invisible(ans) else ans
}
