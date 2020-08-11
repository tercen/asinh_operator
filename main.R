library(tercen)
library(dplyr)

ctx = tercenCtx()
scale             = as.double(ctx$op.value('scale'))
scale_per_channel = as.double(ctx$op.value('scale per channel'))

if (scale_per_channel){
  row_df = ctx %>% rselect() %>% mutate(.ri = row_number()-1)
  scale_name_sym = sym(ctx$rnames[[2]])
  
  ctx %>% 
    select(.ri, .ci, .y) %>% 
    left_join(row_df, by = ".ri") %>% 
    group_by(.ri, .ci) %>% 
    mutate(asinh = asinh(.y)/!!scale_name_sym) %>%
    select(.ri, .ci, asinh) %>%
    ctx$addNamespace() %>%
    ctx$save()
} else {

  ctx  %>%
  select(.ri, .ci, .y) %>%
  group_by(.ri, .ci) %>% 
  mutate(asinh = asinh(.y)/scale) %>%
  select(.ri, .ci, asinh) %>%
  ctx$addNamespace() %>%
  ctx$save()
}