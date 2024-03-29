library(tercen)
library(dplyr)


options("tercen.workflowId"= "237423aa306b89c403ec93aaaa025c4c")
options("tercen.stepId"= "ea0882b5-785f-4254-a715-8c25ce7e4647")
getOption("tercen.workflowId")
getOption("tercen.stepId")



# scale = as.double(ctx$op.value('scale'))
# scale = 5

ctx <- tercenCtx()
scale <- NULL
if (as.character(ctx$op.value('scale')) != "NULL")  scale  = as.integer(ctx$op.value('scale'))


if (is.null(scale)){
  if (length(ctx$rnames) < 2) stop("require to have a scaling value after channel name in projection")
  row_df = ctx %>% rselect() %>% mutate(.ri = row_number()-1)
  scale_colname = sym(ctx$rnames[[2]])
  
  ctx %>% 
    select(.ri, .ci, .y) %>% 
    left_join(row_df, by = ".ri") %>% 
    group_by(.ri, .ci) %>% 
    mutate(asinh = asinh(.y/!!scale_colname)) %>%
    select(.ri, .ci, asinh) %>%
    ctx$addNamespace() %>%
    ctx$save()
} else {
  
  ctx  %>%
    select(.ri, .ci, .y) %>%
    group_by(.ri, .ci) %>% 
    mutate(asinh = asinh(.y/scale)) %>%
    select(.ri, .ci, asinh) %>%
    ctx$addNamespace() %>%
    ctx$save()
}
