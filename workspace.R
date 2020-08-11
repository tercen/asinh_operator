library(tercen)
library(dplyr)

options("tercen.serviceUri"= "https://tercen.com/api/v1/")
options("tercen.username"= "fnaji")
options("tercen.password"= "obama")

options("tercen.workflowId"= "7efcf43fe12c5104b13ccbd7f659226e")
options("tercen.stepId"= "b5abd33b-d136-44e3-9e55-39370ae9b851")
getOption("tercen.workflowId")
getOption("tercen.stepId")



# scale = as.double(ctx$op.value('scale'))
# scale = 5

ctx = tercenCtx()


if (length(ctx$rnames) < 2) stop("require to have a scaling value after channel name in projection")
sn = sym(ctx$rnames[[2]])

row_df = ctx %>% rselect() %>% mutate(.ri = row_number()-1)

(ctx = tercenCtx()) %>% 
  select(.ri, .ci, .y) %>% 
  left_join(row_df, by = ".ri") %>% 
  group_by(.ri, .ci) %>% 
  mutate(asinh = asinh(.y)/!!sn) %>%
  select(.ri, .ci, asinh) %>%
  ctx$addNamespace() %>%
  ctx$save()
  

df <- ctx  %>%
  select(.ri, .ci, .y) %>%
  group_by(.ri, .ci) %>% 
  mutate(asinh = asinh(.y)/scale) %>%
  select(.ri, .ci, asinh) %>%
  ctx$addNamespace() %>%
  ctx$save()
