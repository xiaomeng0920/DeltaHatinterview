#task: delta hat interview questions
#time: 27th April 2025
#name: Xiaomeng Liao


#load packages
pacman::p_load("dplyr","haven", "ggplot2","png","gtsummary","gt","survival","survminer","muhaz","survRM2")


#1. Read data --------------------------------------------------------------------------------------------
data_path<- "C:/Users/XLiao/Downloads/DH task"
setwd("C:/Users/XLiao/Downloads/DH task")
DATA<- read_xpt(file.path(data_path, "adam_dat.xpt"))

data<- DATA %>%
  filter(TRT01P %in%c("tablemab x 52 weeks","vismab x 52 weeks")) %>%
  mutate(event = 1 - CNSR,
         time = AVAL*12)


#2. Patient characteristics 
table(is.na(data$STR01L))
table(is.na(data$STR02L))

#save as pdf
patient_table_pdf<- data %>%
  select(AGE,STR01L,STR02L,TRT01P) %>%
  tbl_summary(by = TRT01P, missing = "no",
              statistic = list(
                AGE ~ "{median}"),
              label = list(AGE ~ "Median age at randomization, years")) %>%
  add_overall()%>%
  as_gt()

gtsave(patient_table_pdf, "Patient characteristics.pdf")

#or save as word document
patient_table_word<- data %>%
  select(AGE,STR01L,STR02L,TRT01P) %>%
  tbl_summary(by = TRT01P, missing = "no",
              statistic = list(
                AGE ~ "{median}"),
              label = list(AGE ~ "Median age at randomization, years")) %>%
  add_overall()%>%
  as_flex_table() 

patient_table_word

doc <- read_docx()
doc <- body_add_flextable(doc, value = patient_table_word)
print(doc, target = "Patient characteristics.docx")


#3. KM curve
range(data$time)
png('KM curve (overall and by treatment).png', res = 300, width = 8, height = 6, units = "in")

#legend
tablemab<- paste("tablemab x 52 weeks","( N=", length(data[data$TRT01P == "tablemab x 52 weeks",]$SUBJID), ")", sep= " ")
vismab<- paste("vismab x 52 weeks","( N=", length(data[data$TRT01P == "vismab x 52 weeks",]$SUBJID), ")", sep= " ")
all_treatment <- paste("overall cohort", "( N=", length(data$SUBJID), ")", sep= " ")

#plot
plot(seq(0,66,1), seq(0,1,1/66) , type="n", xlab="Months", ylab="Survival", axes=F)
axis(2, at=seq(0,1,0.1)) # 2 =left 
axis(1, at=seq(0,66,6))  # 1 = below
grid(nx = NA, ny= NULL)
abline(v = seq(0,66, by=6), lwd = 2, lty = 3, col = "lightgray")
title(main="PFS (overall population and by treatment)", col.main="black", font.main=2)
legend(50,1, c(tablemab, vismab, all_treatment), cex= 0.55, col=c("#9C0824","#1C63A1", "#3A3A3A"), lwd=c(2,2,2), lty=1)

#fit
data_T<- data[data$TRT01P == "tablemab x 52 weeks",]
data_v<- data[data$TRT01P == "vismab x 52 weeks",]

fit_tablemab<- survfit(Surv(data_T$time,data_T$event) ~ 1)
fit_vismab<- survfit(Surv(data_v$time,data_v$event) ~ 1)
fit<- survfit(Surv(data$time,data$event) ~ 1)

lines(fit_tablemab,  mark.time =T, lty=1, col="#9C0824", conf.int=F, lwd=1)
lines(fit_vismab,  mark.time =T, lty=1, col="#1C63A1", conf.int=F, lwd=1)
lines(fit,  mark.time =T, lty=1, col="#3A3A3A", conf.int=F, lwd=1)
dev.off()

#median survival
medians_tablemab <- summary(fit_tablemab)$table["median"]
medians_tablemab

medians_vismab <- summary(fit_vismab)$table["median"]
medians_vismab

medians_all <- summary(fit)$table["median"]
medians_all

#question 3c
data_3c<- data %>%
  filter(TRT01P == "tablemab x 52 weeks")%>%
  filter(STR01 == "Negative")

fit_3c<- survfit(Surv(data_3c$time,data_3c$event) ~ 1)
medians_3c<- summary(fit_3c)$table["median"]
medians_3c

#4.Hazard
png('Smooth hazard plot.png', res = 300, width = 10, height = 6, units = "in")

haz_T <- muhaz(data_T$time, data_T$event)
haz_v <- muhaz(data_v$time, data_v$event)

plot(haz_T, col = "#AE123A", lwd = 2, xlab = "Time", ylab = "Hazard", main = "Smoothed Hazard plot")
lines(haz_v, col = "#26456E", lwd = 2)

legend("bottomright", legend = c("tablemab x 52 weeks", "vismab x 52 weeks"), col = c("#AE123A", "#26456E"), lwd = 2,cex= 0.7)
dev.off()

#smooth hazard plot shows the instantaneous event risk versus time. Both treatments showed similar trends in terms of hazard: the risk 
#was very low initially and then went up quickly until around 20 weeks; from 20 to 40 weeks, the increasing speed slows down
#,the hazard increase slowly and remains in a more steady level. After 45 weeks, the hazard increase quickly again.
#even though the event risk is higher for vismab in the first 5 weeks, during most of the time (5 to 45 weeks), patients receiving tablemab experienced a slightly higher risk of event than patients in vismab group.
#the difference between them becomes smaller after week 40 and disappear after 45 weeks as the two curves crossed.
#proportional hazard assumption might not hold as the hazard difference constantly changed over time and the two curves crossed twice.

#5.	Cox model
cox_model <- coxph(Surv(time, event) ~ TRT01P, data = data)
summary(cox_model)

check_ph<- cox.zph(cox_model)
check_ph
plot(check_ph)

#p value for treatment and the whole model is 0.82, which is much larger than 0.05;
#looking at the the Schoenfeld residuals plot, the scatter point do not show any trend, the solid line is almost horizontal
#These all indicate the proportional hazard assumption hold, thus, the use of Cox regression model is reasonable.
#from the model, hazard ratio between vismab and tablemab is 0.91,95% CI (0.75,1.10)
#the patients in vismab group had a 9% lower risk of disease progression or death compared to patients receiving tablemab, however, 
#this difference is not statistically significant, there is 95% probability that vismab patients had 25% lower risk ir 10% higher chance of disease progression or death.

#6.the restricted mean survival time
tau <- 48
data$treatment_numeric <- ifelse(data$TRT01P == "tablemab x 52 weeks", 1, 0)

rmst_result <- rmst2(time = data$time, status = data$event, arm = data$treatment_numeric, tau = 48)
rmst_result

#restricted mean survival gives the average survival time until particular time.
#at 48 months, the restricted mean survival time was 31.5 months for patients in tablemab group and 32.7 months for vismab patients.
#there is a non-statistically significant difference of -1.2 months (95% CI: -3.4 to 1.0; p=0.298) between the groups, indicating
#patients receiving tablemab had a average of 1.2 months lower survival time than patients in vismab group. This difference is not statistically meaningful.


