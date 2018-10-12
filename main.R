library(tercen)
library(dplyr)

scale = as.double(ctx$op.value('scale'))

(ctx = tercenCtx()) %>% 
  select(.y) %>% 
  transmute(asinh = asinh(.y)/scale) %>%
  ctx$addNamespace() %>%
  ctx$save()