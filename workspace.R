library(tercen)
library(dplyr)

options("tercen.serviceUri" = "http://127.0.0.1:5400/api/v1")
options("tercen.username"= "admin")
options("tercen.password" = "admin")

options("tercen.workflowId"= "1b5b9c21c7cbf0613462606fba00b2f2")
options("tercen.stepId"= "3-1")
getOption("tercen.workflowId")
getOption("tercen.stepId")


# scale = as.double(ctx$op.value('scale'))
scale = 5

(ctx = tercenCtx()) %>% 
  select(.y) %>% 
  transmute(asinh = asinh(.y)/scale) %>%
  ctx$addNamespace() %>%
  ctx$save()
