###------------------------------------------------------------------------
###                         MY CUSTOM FUNCTIONS                         ---
###                          MICHAL FOLWARCZNY                          ---
###                           MICHALF@TUTA.IO                           ---
###------------------------------------------------------------------------

###------------------------------------------------------------------------
###                              COHEN'S D                              ---
###------------------------------------------------------------------------

coh.d <- function(x, y, ci = .95, digits = 2) {
  lx <- length(x) - 1
  ly <- length(y) - 1
  diff <- abs(mean(x, na.rm = T) - mean(y, na.rm = T))
  d1 <- lx * var(x, na.rm = T) + ly * var(y, na.rm = T)
  d2 <- d1 / (lx + ly)
  d3 <- sqrt(d2)
  d <- diff / d3
  lower.ci <- psych::cohen.d.ci(d,
                                n1 = length(x),
                                n2 = length(y),
                                alpha = 1-ci)[1]
  upper.ci <- psych::cohen.d.ci(d,
                                n1 = length(x),
                                n2 = length(y),
                                alpha = 1-ci)[3]
  print(round(data.frame("Cohen d" = d,
                         "CI.lower" = lower.ci,
                         "CI.Upper" = upper.ci,
                         "CI.level" = ci),
              digits))
}

# Example: coh.d(df1$sum[df1$group == "control"],
# df1$sum[df1$group == "experimental"], ci = .95, digits = 3)

###-------------------------------------------------------------------------
###                             DESCRIPTIVES                             ---
###-------------------------------------------------------------------------
ds <- function(x, ci = .95, digits = 2) {
  print(round(data.frame("Mean" = mean(x),
                         "SD" = sd(x),
                         "CI.Lower" = mean(x) -
                           qt(ci + (1 - ci)/2,
                              df = length(x) - 1) * sd(x) / sqrt(length(x)),
                         "CI.Upper" = mean(x) +
                           qt(ci + (1 - ci)/2,
                              df = length(x) - 1) * sd(x) / sqrt(length(x)),
                         "SE" = sd(x, na.rm = T)/sqrt(length(x[!is.na(x)])),
                         "IQR" = IQR(x),
                         "Q1" = quantile(x, 0.25)[[1]],
                         "Q3" = quantile(x, 0.75)[[1]],
                         "Min" = min(x),
                         "Max" = max(x),
                         "Median" = median(x),
                         "CI.level" = ci), digits))
  par(mfrow = c(2,2), cex = .6)
  hist(x, col = "lightsteelblue1")
  abline(v = mean(x),
         col = "royalblue",
         lwd = 2,
         lty = 2)
  abline(v = median(x),
         col = "chocolate3",
         lwd = 2,
         lty = 2)
  legend(x = "topright",
         c("MEAN", "MEDIAN"),
         col = c("chocolate3", "royalblue"),
         lwd = c(5, 5),
         bty = "n")
  dens <- density(x)
  plot(dens, lwd = 1.5)
  polygon(dens, col = "lightsteelblue1")
  abline(v = mean(x),
         col = "royalblue",
         lwd = 2,
         lty = 2)
  abline(v = median(x),
         col = "chocolate3",
         lwd = 2,
         lty = 2)
  legend(x = "topright",
         c("MEAN", "MEDIAN"),
         col = c("chocolate3", "royalblue"),
         lwd = c(5, 5),
         bty =  "n")
  boxplot(x, col = "lightsteelblue1", horizontal = T)
  qqnorm(x)
  qqline(x, col = "powderblue", lwd = 2)
  par(mfrow = c(1,1))
}
# Example: ds(df1$age)

###-------------------------------------------------------------------------
###                    DOWNLOADING & LOADING PACKAGES                    ---
###-------------------------------------------------------------------------
ipak <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if(length(new.pkg)) 
    install.packages(new.pkg, dependencies = T)
  sapply(pkg, require, character.only = T) 
}

###-------------------------------------------------------------------------
###                               Z-VALUES                               ---
###-------------------------------------------------------------------------
z.val <- function(score, dataset) {
  z.val <- (score - mean(dataset)) / sd(dataset)
  round(z.val, digits = 2) }



###-------------------------------------------------------------------------
###                               Benford's Law                          ---
###-------------------------------------------------------------------------
# Benford's Law for checking, whether data distributions follows that found in real life

# benlaw <- function(d) log10(1 + 1 / d)
# digits <- 1:9
# 
# jpeg(file = "test1.jpeg")
# 
# baseBarplot <- barplot(benlaw(digits), names.arg = digits, xlab = "First Digit", 
#                        ylim = c(0, .35))
# firstDigit <- function(x) substr(gsub('[0.]', '', x), 1, 1)
# pctFirstDigit <- function(x) data.frame(table(firstDigit(x)) / length(x))
# 
# test1 <- pctFirstDigit(df1$column1)
# 
# lines(x = baseBarplot[,1], y = test1$Freq, col = "red", lwd = 4, 
#       type = "b", pch = 23, cex = 1.5, bg = "red")
# 
# dev.off()

###-------------------------------------------------------------------------
###                           OUTLIERS REMOVAL                           ---
###-------------------------------------------------------------------------
# remove.outliers <- function(x, na.rm = TRUE, ...) {
#   qnt <- quantile(x, probs=c(.25, .75), na.rm = na.rm, ...)
#   H <- 1.5 * IQR(x, na.rm = na.rm)
#   y <- x
#   y[x < (qnt[1] - H)] <- NA
#   y[x > (qnt[2] + H)] <- NA
#   y
# }
# example
# df1.outrem1 <- df1
# df1.outrem1$lux.mean <- remove_outliers(df1$lux.mean)
# df1.outrem1$attitude.mean <- remove_outliers(df1$attitude.mean)
# df1.outrem1 <- drop_na(df1.outrem1)
# df1.outrem1

###------------------------------------------------------------------------
###                       LOOP FOR LM ACCROSS DVS                       ---
###------------------------------------------------------------------------
# # Select column (DV) names
# vars <- names(df1)[4:11]
# # Prepare an empty container for results
# res <- rep(NA, length(vars))
# # Write a loop over columns
# for(i in 1:length(vars)) {
#   # Write the formula, left side = DV name, right side = IV name
#   myformula <- paste(vars[i],'~ hair.color')
#   # Fit the linear model
#   mod <- lm(myformula, data = df1)
#   # Save results in an empty container; here, we extract p-values
#   res[i] <- summary(mod)$coefficients[2,4]
# }
# # Print the results
# res

###------------------------------------------------------------------------
###                   EFFECT SIZES IN CHI-SQUARE TEST                   ---
###------------------------------------------------------------------------
# cramer.v <- function(x, y, n) {
#   phi <- sqrt(as.numeric(chisq.test(table(x, y))$statistic)/
#                 n * (min(dim(table(x, y))) -1))
#   round((phi), digits = 3) 
# }
# Example: phi.chisq(df1$car, df1$group, 101)

###------------------------------------------------------------------------
###                       RUN BEFORE EACH SESSION                       ---
###------------------------------------------------------------------------

# Banners for commenting code chunks
# library(bannerCommenter)
# banner("Run before each session",
#        emph = T, fold = T, center = T, bandChar = "-", snug = F,
#        numLines = 1, leftSideHashes = 3, rightSideHashes = 3, maxChar = 80)

# Load custom functions
# source(file = "C:/Users/fmich/Desktop/Sync/.Misc/R_misc/MyFuncs.R")

# Set working directory
# setwd(dirname(rstudioapi::getActiveDocumentContext()$path))