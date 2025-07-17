library(psych)
library(polycor)
library(snc)


dassfull<-read.csv("E:/SNC/dass21_subset.csv")

R <- hetcor(dassfull)$correlations

psych::alpha(R)

# Subscale item names
dep_items <- c("Q3A", "Q13A", "Q16A", "Q24A", "Q25A", "Q28A", "Q30A")
anx_items <- c("Q2A", "Q7A", "Q15A", "Q19A", "Q20A", "Q23A", "Q27A")
str_items <- c("Q1A", "Q6A", "Q8A", "Q11A", "Q12A", "Q14A", "Q18A")

# Alpha per subscale
psych::alpha(dassfull[dep_items])
psych::alpha(dassfull[anx_items])
psych::alpha(dassfull[str_items])

# SNC
R_dep <- hetcor(dassfull[dep_items])$correlations
R_anx <- hetcor(dassfull[anx_items])$correlations
R_str <- hetcor(dassfull[str_items])$correlations

snc(R_dep)
snc(R_anx)
snc(R_str)

#psych::alpha(dassfull[str_items],check.keys=TRUE)

# Reorder columns for factor alignment
ordered_items <- c(dep_items, anx_items, str_items)
dass_ordered <- dassfull[ , ordered_items]

# Create vector of factor labels
factors <- c(rep("Depression", 7), rep("Anxiety", 7), rep("Stress", 7))

# Run SNC with factor labeling
snc_dass <- snc(R, k = 2, factors = factors)

# Print the results
print(snc_dass)


####################
# Roesenberg Scale #
####################

rosfull<-ros <- read.delim("E:/SNC/data.csv", sep = "\t", stringsAsFactors = FALSE)

head(rosfull)
ros<-rosfull[,1:10]
head(ros)

ros2 <- hetcor(ros)$correlations

snc(ros2, k = 2)

psych::alpha(ros2)

alpha(ros2,check.keys=TRUE)