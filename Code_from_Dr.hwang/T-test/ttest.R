library(dplyr)

# The data set
d <- read.csv("NIBR-H-19_pH, salinity test_211215.csv")
#dt <- as.data.frame(t(as.matrix(d)))
myData <- as.data.frame(d) %>% select(LB_ph, pH2) #pH: pH3, pH4 ...
#myData <- as.data.frame(d) %>% select(LB_per, Percent0) #salinity: Percent1, Percent2 ...

# Step 1. independent? : yes
# Step 2. follows a normal distribution? : test by "Shapiro" ; p > 0.05 = normal (a parametric alternative) / p < 0.05 = not normal (a non-parametric alternative: Wilcoxon rank sum test or Mann-Whitney test)
with(myData, shapiro.test(LB_ph))
# Step 3. same variances? : test by "var.test (F-test)" ; p > 0.05 = no different (var.equal = TRUE)  / p < 0.05 = different (var.equal = FALSE)
var.test(myData$LB_ph, myData$pH2, data = myData)
# Step 4. t.test: p < 0.05 = different the two group
t.test(myData$LB_ph, myData$pH2, alternative="two.sided", var.equal = TRUE)











##### ------------------------------ Reference ---------------------------------
# Create a data frame
T1 <- c(0.02, 0.018) #37oC
T2 <- c(0.054, 0.064) #33
T3 <- c(0.079, 0.111) #30

myData <- data.frame(
  group = rep(c("T1", "T3"), each = 2),
  weight = c(T1, T3))

print(myData)



##### 1) unpaired two-samples t-test: a parametric alternative (test whether if average is 0)
# http://www.sthda.com/english/wiki/unpaired-two-samples-t-test-in-r

# Step 1. independent? : yes

# Step 2. follows a normal distribution? : test by "Shapiro" ; p > 0.05 = normal / p < 0.05 = not normal
with(myData, shapiro.test(weight[group == "T1"])) #sample size must be between 3 and 5000
with(myData, shapiro.test(weight[group == "T2"]))

# Step 3. same variances? : test by "var.test (F-test)" ; p > 0.05 = no different (var.equal = TRUE)  / p < 0.05 = different (var.equal = FALSE)
var.test(weight ~ group, data = myData)

# Step 4. t.test: p < 0.05 = different the two group
t.test(weight ~ group, data = myData, alternative="two.sided", var.equal = TRUE)

# 37도 vs 33도: p-value = 0.01586
# 37도 vs 30도: p-value = 0.04173


##### 2) unpaired two-samples t-test: a non-parametric alternative (test whether if median is 0)
# http://www.sthda.com/english/wiki/unpaired-two-samples-wilcoxon-test-in-r

# Wilcoxon rank sum test or Mann-Whitney test: p > 0.05 = no different / p < 0.05 = different
wilcox.test(T1, T2) #p=0.3333
wilcox.test(T1, T3) #p=0.3333
wilcox.test(T2, T3) #p=0.3333


