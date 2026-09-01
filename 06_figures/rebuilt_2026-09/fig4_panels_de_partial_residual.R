suppressMessages({library(data.table); library(ggplot2); library(splines)})
setwd("/Volumes/S840/04mypaper/lin/conflictresolution/fusio2trait")
d <- fread("model_input_1678.tsv")
COLS <- c(Male="#2C5985", Female="#C0392B"); GREY <- "#6E7B85"
TH <- theme_classic(base_size=10)+theme(axis.text=element_text(size=8),
      axis.title=element_text(size=9), plot.subtitle=element_text(size=8))
FULL <- lp ~ factor(chap)+rho_GE_male+rho_GE_female+lmort+ns(onset_all,3)
mf <- lm(FULL, d)
pr <- function(drop){ m0 <- update(mf, as.formula(paste(".~.-", drop)))
  (summary(mf)$r.squared-summary(m0)$r.squared)/(1-summary(m0)$r.squared) }

X <- model.matrix(mf); cf <- coef(mf)
sc <- grep("onset_all", colnames(X), fixed=TRUE)
comp <- as.vector(X[, sc, drop=FALSE] %*% cf[sc])
pd <- data.table(onset=d$onset_all, y=comp - mean(comp) + resid(mf))
pD <- ggplot(pd, aes(onset, y)) +
  geom_point(colour=GREY, size=0.7, alpha=0.28, stroke=0) +
  geom_smooth(method="lm", formula=y~ns(x,3), se=FALSE, colour=GREY, linewidth=0.9) +
  labs(x="Median age at first recorded event (years)", y="Prevalence (partial residual)",
       subtitle=sprintf("partial R² = %.3f", pr("ns(onset_all, 3)"))) + TH

OTH <- "factor(chap) + lmort + ns(onset_all,3)"
avdt <- function(v,lab){ o <- setdiff(c("rho_GE_male","rho_GE_female"), v)
  f0 <- paste("~", OTH, "+", o)
  data.table(rx=resid(lm(as.formula(paste(v,f0)),d)), ry=resid(lm(as.formula(paste("lp",f0)),d)), axis=lab)}
E <- rbind(avdt("rho_GE_male","Male"), avdt("rho_GE_female","Female"))
joint <- pr("rho_GE_male-rho_GE_female")
sE <- data.table(axis=c("Male","Female"), pp=c(pr("rho_GE_male"), pr("rho_GE_female")))
sE[, lab := sprintf("'%s: partial '*italic(R)^2=='%.3f'", substr(axis,1,1), pp)]
jlab <- sprintf("'joint: partial '*italic(R)^2=='%.3f'", joint)
cat(sprintf("joint block partial R2 = %.4f\n", joint))
print(sE)
pE <- ggplot(E, aes(rx, ry, colour=axis)) +
  geom_point(size=0.7, alpha=0.30, stroke=0) +
  geom_smooth(aes(linetype=axis), method="lm", formula=y~x, se=FALSE, linewidth=0.8) +
  geom_text(data=sE[axis=="Male"],   aes(x=-Inf,y=Inf,label=lab), hjust=-0.05, vjust=1.5, size=2.4, parse=TRUE, show.legend=FALSE, colour=COLS["Male"]) +
  geom_text(data=sE[axis=="Female"], aes(x=-Inf,y=Inf,label=lab), hjust=-0.05, vjust=3,   size=2.4, parse=TRUE, show.legend=FALSE, colour=COLS["Female"]) +
  annotate("text", x=-Inf, y=Inf, label=jlab, hjust=-0.05, vjust=4.5, size=2.4, parse=TRUE, colour="grey25") +
  scale_colour_manual(values=COLS, name=NULL, labels=c(Male="Male fertility", Female="Female fertility")) +
  scale_linetype_manual(values=c(Male="solid", Female="22"), guide="none") +
  labs(x="Disease–fertility coupling (residual)", y="Prevalence (residual)") +
  guides(colour=guide_legend(override.aes=list(size=2, alpha=1))) + TH +
  theme(legend.position="inside", legend.position.inside=c(0.80,0.14),
        legend.text=element_text(size=7.5), legend.key.size=unit(3.5,"mm"), legend.background=element_blank())
dev <- if (capabilities("cairo")) cairo_pdf else pdf
ggsave("figures_final/Fig4d_onset_av.pdf",    pD, width=89, height=80, units="mm", device=dev, bg="white")
ggsave("figures_final/Fig4e_coupling_av.pdf", pE, width=89, height=80, units="mm", device=dev, bg="white")
cat(sprintf("onset partial R2 = %.4f\n", pr("ns(onset_all, 3)")))
