library(devtools)
document("ZINB")
load_all("ZINB")

build("ZINB")
check("ZINB", manual = T)