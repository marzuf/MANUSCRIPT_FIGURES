config_dir ="."

#' Read params from files
init_params = function(fn, prev_params=NULL)
{
    t = read.table(fn, sep = "=", fill = T, strip.white=T, stringsAsFactors=FALSE, quote="")
    params = as.character(t[, 2])
    names(params) = t[, 1]

    if (any(duplicated(names(params))))
    {
        warning("duplicated params in conf file! check ", fn)
        return(NA)
    }

    if (!is.null(prev_params)){
        if(any(names(prev_params)) %in% names(params)){
		stop(sprintf("duplicaterd parameters: %s", names(prev_parames)[names(prev_params) %in% names(params)]))
        }
    }
    params <- as.list(params)
    # interpret some values as R exprssions
    expr_idx <- grep("\\@R$", names(params))
    # Strip names from R
    params[expr_idx] <- lapply(params[expr_idx], function(x) eval(parse(text=x)))
    names(params)[expr_idx] <- sub("\\@R$", "", names(params[expr_idx]))

    return(params)
}


#' Dump config file templates.
#'
#' Dump templates of config files required by the package.
#'
#' @param config_dir Directory to dump to files to.
#'
#'

.shaman_dump_config <- function(config_dir)
{
    config_files <- dir(system.file("config", package='shaman'), full.names=T)
    dir.create(config_dir, recursive=T, showWarnings=FALSE)
    ret <- file.copy(config_files, config_dir, recursive=T, overwrite=FALSE)
    if (!all(ret))
    {
        warning("couldn't dump config files to ", config_dir, "\n  Perhaps they're already there? ")
    } else
    {
        message("Dumped config files to: ", config_dir)
    }
}

#' Load all the configuration files required for the package.
#'
#' Loads the configuration files which were filled by the user. Accepts either directory path
#' from \code{config_dir}, or individual files supplied by the other parameters.
#'
#' @param config_dir A directory with all the defined config files. Should have the same
#' files as those that were exported with \code{shaman_dump_config}.
#' @param shaman_config Parameter file
#'
#'
.shaman_load_config <- function(config_dir, shaman_config=file.path(".", "shaman.conf"), reset=FALSE)
{
    op <- options()
    all_paths <- list("shaman.shaman_conf_path"=shaman_config)

    op.shaman <- init_params(shaman_config)
    op.shaman<- c(op.shaman, all_paths)

    # Check if some params were saved already
    already_set <- names(op) %in% names(op.shaman)
    if (any(already_set) & reset==TRUE)
    {
        message(sprintf("%i parameters were already defined in this sessions and they will be overwritten!",
                sum(already_set)))
    }

    if (!reset) {
	op.shaman = op.shaman[!already_set]
    }
    options(op.shaman)

#    gset_input_mode(autocompletion = FALSE, interactive = FALSE)
    invisible()
}


#######################################################################
# .shaman_check_config
#
# internal function for checking that the configuration has been loaded
########################################################################
.shaman_check_config <- function(params) {
  if (!all( params %in% names(options()) ))
    {
        .shaman_load_config(system.file("config", package="shaman"))
    }
}
