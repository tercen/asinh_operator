library(tercenApi)
library(tercen)
library(data.table)
library(dtplyr)
library(dplyr, warn.conflicts = FALSE)

ctx <- tercenCtx()
scale = ctx$op.value("scale", type = as.integer, default = 5)

if (scale < 0) {
  if (length(ctx$rnames) < 2)
    stop("require to have a scaling value after channel name in projection")
  
  row_df =
    ctx$rselect(ctx$rnames[[2]]) %>% mutate(.ri = row_number() - 1) %>% lazy_dt()
  
  scale_colname <- sym(ctx$rnames[[2]])
  
  ctx %>%
    select(.ri, .ci, .y) %>%
    lazy_dt() %>%
    left_join(row_df, by = ".ri") %>%
    mutate(asinh = asinh(.y / !!scale_colname)) %>%
    select(.ri, .ci, asinh) %>%
    as_tibble() %>%
    ctx$addNamespace() %>%
    ctx$save()
} else {
  ctx  %>%
    select(.ri, .ci, .y) %>%
    lazy_dt() %>%
    mutate(asinh = asinh(.y / scale)) %>%
    select(.ri, .ci, asinh) %>%
    as_tibble() %>%
    ctx$addNamespace() %>%
    ctx$save()
}
