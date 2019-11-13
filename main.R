library(tercen)
library(dplyr)

ctx = tercenCtx()
scale = as.double(ctx$op.value('scale'))

ctx %>% 
  select(.y) %>% 
  transmute(asinh = asinh(.y/scale)) %>%
  ctx$addNamespace() %>%
  ctx$save()