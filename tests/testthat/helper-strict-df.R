.strict_df <- function(x) {
  class(x) <- c("strict_df", class(x))
  x
}

# The `$` of a tibble, which warns on an absent column, without tibble.
registerS3method("$", "strict_df", function(x, name) {
  if (!name %in% names(x)) {
    warning(sprintf("Unknown or uninitialised column: `%s`.", name),
            call. = FALSE)
  }
  .subset2(x, name)
})
