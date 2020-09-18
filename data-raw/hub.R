# download hub resources locally
library(AnnotationHub)
AnnotationHub(cache = here::here("data"))[["AH73905"]]
AnnotationHub(cache = here::here("data"))[["AH33458"]]
