#### library ####
library('psych')
library('irr')
library("dplyr")
library("tidyverse")
library("readxl")
library("corrplot")
library("RColorBrewser")
library("rstatix")
library("ggpubr")
library("ggbeeswarm")
library("openxlsx")
library("writexl")
library("epiR")
library('irrCAC')


#### setwd ####
setwd("...")

pt_demography = read_excel("...",
                            sheet = "Pt")

pt_anxiety_survey = read_excel("...",
                               sheet = "Anxiety")

TB_sup_vs_pro = read_excel ("...", sheet = "PI-QUAL",
                            range = "A3:Y55")

IC_sup_vs_pro = read_excel ("...", sheet = "PI-QUAL",
                            range = "A3:Y48")

colnames(TB_sup_vs_pro)

TB_rectal_air = read_excel("...", sheet = "air",
                           range = "A1:E54")

Buscopan_data = read_excel ("...")

supine_air_quan = read_excel("...")



#### process data ####
pt_demography = convert_dates(pt_demography, "ScanDate")

pt_demography_prone = pt_demography %>%
  filter(is.na(Prone))


pt_anxiety_survey = convert_dates(pt_anxiety_survey, "ScanDate")

pt_anxiety_survey_prone = 
  pt_anxiety_survey %>%
  filter (ExamID %in% pt_demography_prone$ExamID)

pt_anxiety_survey_prone =
  pt_anxiety_survey_prone %>%
  mutate(Q16 = as.numeric(Q16))

# Copy values from "Sup.T2.Sag" to "Pro.T2.Sag"
IC_sup_vs_pro$Pro.T2.Sag <- IC_sup_vs_pro$Sup.T2.Sag

#### process data with predefined function in "01-01_Supine_vs_prone_PIQUAL_20241024.R" ####
IC_sup_vs_pro = convert_dates(IC_sup_vs_pro, "ScanDate")

IC_sup_vs_pro = processData_sum_individual_seq(IC_sup_vs_pro)

IC_sup_vs_pro = processData_PIQUALv.2_score(IC_sup_vs_pro)

IC_sup_vs_pro = processData_PIQUALv.2_individual_seq_to_01(IC_sup_vs_pro)

####### part 1 #######
#### pt demographics summary ####
summary(pt_demography_prone)

pt_demography_prone %>%
  get_summary_stats(Age, type = "quantile")

pt_demography_prone %>%
  get_summary_stats(PSA, type = "quantile")

pt_demography_prone %>%
  get_summary_stats(PSAd, type = "quantile")


pt_demography_prone %>%
  count(`PI-RADS`, name = "Frequency") %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)

table(pt_demography_prone$`PI-RADS`) 

pt_demography_prone %>%
  count(Patho, name = "Frequency")


pt_anxiety_survey_prone %>%
  count(Q16, name = "Frequency") %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)


pt_anxiety_survey_prone %>%
  filter(!is.na(Q16)) %>%
  summarise(non_na_count = n())

pt_anxiety_survey_prone %>%
  get_summary_stats(Q16, type = "quantile")


#### TB_rectal_air supine vs prone summary (have to run 01-01 first) ####
TB_rectal_air %>%
  count(Supine.Rectal.Score, name = "Frequency") %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)

TB_rectal_air %>%
  mutate(Group = case_when(
    Supine.Rectal.Score %in% 1:3 ~ "Group 1-3",
    Supine.Rectal.Score %in% 4:5 ~ "Group 4-5"
  )) %>%
  count(Group, name = "Frequency") %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)


TB_rectal_air %>%
  count(Prone.Rectal.Score, name = "Frequency") %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)

TB_rectal_air %>%
  mutate(Group = case_when(
    Prone.Rectal.Score %in% 1:3 ~ "Group 1-3",
    Prone.Rectal.Score %in% 4:5 ~ "Group 4-5"
  )) %>%
  count(Group, name = "Frequency") %>%
  mutate(Percentage = Frequency / sum(Frequency) * 100)


#### supine vs prone T2, DWI summary (have to run 01-01 first) ####
TB_sup_vs_pro %>%
  select(Sup.T2_total, Pro.T2_total, Sup.DWI_total, Pro.DWI_total) %>%
  get_summary_stats(type = "full")  # "full" gives all available stats

####### part 2 #######

# Copy values from "Sup.T2.Sag" to "Pro.T2.Sag"
TB_sup_vs_pro$Pro.T2.Sag <- TB_sup_vs_pro$Sup.T2.Sag

# Convert specified columns to numeric
cols_to_convert <- c("Sup.T2.Artef", "Sup.T2.SNR", "Sup.T2.Struct", "Sup.T2.Sag", 
                     "Sup.DWI.b_SNR", "Sup.DWI.ADC", "Sup.DWI.Artef", "Sup.DWI.5mm", 
                     "Sup.DCE.Enh", "Sup.DCE.Anat", "mp.Sup.Total_10", "mp.Sup.Overall_3", 
                     "Pro.T2.Artef", "Pro.T2.SNR", "Pro.T2.Struct", "Pro.T2.Sag", 
                     "Pro.DWI.b_SNR", "Pro.DWI.ADC", "Pro.DWI.Artef", "Pro.DWI.5mm", 
                     "mp.Pro.Total_10", "mp.Pro.Overall_3")

TB_sup_vs_pro[cols_to_convert] <- lapply(TB_sup_vs_pro[cols_to_convert], as.numeric)

TB_sup_vs_pro %>% 
  select (ScanDate) %>%
  print(n=30)


# change scan date to the correct format
TB_sup_vs_pro <- TB_sup_vs_pro %>%
  mutate(ScanDate = case_when(
    !is.na(as.numeric(ScanDate)) ~ as.Date(as.numeric(ScanDate), origin = "1899-12-30"),
    grepl("^\\d{2}/\\d{2}/\\d{4}$", ScanDate) ~ as.Date(ScanDate, format = "%d/%m/%Y"),
    TRUE ~ as.Date(NA)  # Assign NA for any values that do not match the expected patterns
  ))


TB_sup_vs_pro %>% 
  select (ScanDate) %>%
  print(n=30)

Buscopan_data <- Buscopan_data %>%
  mutate(ScanDate = case_when(
    !is.na(as.numeric(ScanDate)) ~ as.Date(as.numeric(ScanDate), origin = "1899-12-30"),
    grepl("^\\d{2}/\\d{2}/\\d{4}$", ScanDate) ~ as.Date(ScanDate, format = "%d/%m/%Y"),
    TRUE ~ as.Date(NA)  # Assign NA for any values that do not match the expected patterns
  ))



#### Handy function ####

convert_dates <- function(data, date_column, origin_date = "1899-12-30") {
  
  data %>%
    mutate(!!sym(date_column) := case_when(
      !is.na(as.numeric(!!sym(date_column))) ~ as.Date(as.numeric(!!sym(date_column)), origin = origin_date),
      grepl("^\\d{2}/\\d{2}/\\d{4}$", !!sym(date_column)) ~ as.Date(!!sym(date_column), format = "%d/%m/%Y"),
      TRUE ~ as.Date(NA)  # Assign NA for values not matching patterns
    ))
}

processData_sum_individual_seq = function(data) {
  result = data %>%
    mutate(Sup.T2_total = rowSums(select(., Sup.T2.Artef:Sup.T2.Sag), na.rm = FALSE)) %>%
    mutate(Sup.DWI_total = rowSums(select(., Sup.DWI.b_SNR:Sup.DWI.5mm), na.rm = FALSE)) %>%
    mutate(Sup.DCE_total = rowSums(select(., Sup.DCE.Enh:Sup.DCE.Anat), na.rm = FALSE)) %>%
    mutate(Pro.T2_total = rowSums(select(., Pro.T2.Artef:Pro.T2.Sag), na.rm = FALSE)) %>%
    mutate(Pro.DWI_total = rowSums(select(., Pro.DWI.b_SNR:Pro.DWI.5mm), na.rm = FALSE)) %>%
    select(No, AnonID, ScanDate,
           starts_with("Sup.T2."), Sup.T2_total,
           starts_with("Sup.DWI."), Sup.DWI_total,
           starts_with("Sup.DCE."), Sup.DCE_total,
           mp.Sup.Total_10, mp.Sup.Overall_3,
           starts_with("Pro.T2."), Pro.T2_total,
           starts_with("Pro.DWI."), Pro.DWI_total,
           mp.Pro.Total_10, mp.Pro.Overall_3)
  return(result)
}


processData_PIQUALv.2_score = function(data) {
  data %>%
    mutate(mp.Sup.Overall_3.check = 
             ifelse((Sup.T2_total) == 4 & (Sup.DWI_total) == 4 & (Sup.DCE_total) == 2, 3,
                    ifelse((Sup.T2_total) == 4 & (Sup.DWI_total) == 4 & (Sup.DCE_total) < 2, 2,
                           ifelse((Sup.T2_total) >= 3 & (Sup.DWI_total) >= 3, 2, 
                                  ifelse((Sup.DCE_total) == 2 & ( (Sup.T2_total) == 4 | (Sup.DWI_total) == 4), 2,
                                         1
                                  ))))
    ) %>%
    
    mutate(mp.Pro.Overall_3.check = 
             ifelse((Pro.T2_total) == 4 & (Pro.DWI_total) == 4 & (Sup.DCE_total) == 2, 3,
                    ifelse((Pro.T2_total) == 4 & (Pro.DWI_total) == 4 & (Sup.DCE_total) < 2, 2,
                           ifelse((Pro.T2_total) >= 3 & (Pro.DWI_total) >= 3, 2, 
                                  ifelse((Sup.DCE_total) == 2 & ( (Pro.T2_total) == 4  | (Pro.DWI_total) == 4), 2,
                                         1
                                  ))))
    ) %>%
    
    mutate(mp.Sup.T2.Pro.DWI.Overall_3.check = 
             ifelse((Sup.T2_total) == 4 & (Pro.DWI_total) == 4 & (Sup.DCE_total) == 2, 3,
                    ifelse((Sup.T2_total) == 4 & (Pro.DWI_total) == 4 & (Sup.DCE_total) < 2, 2,
                           ifelse((Sup.T2_total) >= 3 & (Pro.DWI_total) >= 3, 2, 
                                  ifelse((Sup.DCE_total) == 2 & ((Sup.T2_total) == 4  | (Pro.DWI_total) == 4), 2,
                                         1
                                  ))))
    ) %>%
    
    mutate(bp.Sup.Overall_3.check = 
             ifelse((Sup.T2_total) == 4 & (Sup.DWI_total) == 4, 3,
                    ifelse((Sup.T2_total) >= 3 & (Sup.DWI_total) >= 3, 2, 
                           1))
    ) %>%
    
    mutate(bp.Pro.Overall_3.check = 
             ifelse((Pro.T2_total) == 4 & (Pro.DWI_total) == 4, 3,
                    ifelse((Pro.T2_total) >= 3 & (Pro.DWI_total) >= 3, 2, 
                           1))
    ) %>%
    
    mutate(bp.Sup.T2.Pro.DWI.Overall_3.check = 
             ifelse((Sup.T2_total) == 4 & (Pro.DWI_total) == 4, 3,
                    ifelse((Sup.T2_total) >= 3 & (Pro.DWI_total) >= 3, 2, 
                           1))
    )
  
}

### Convert V.2 individual sequence to (0,1) binary system
processData_PIQUALv.2_individual_seq_to_01 = function(data) {
  data %>%
    mutate(Sup.T2_01 = 
             ifelse(Sup.T2_total >= 3, 1, 0)) %>%
    mutate(Sup.DWI_01 = 
             ifelse(Sup.DWI_total >= 3, 1, 0)) %>%
    mutate(Sup.DCE_01 = 
             ifelse(is.na(Sup.DCE_total), NA, 
                    ifelse(Sup.DCE_total == 2, 1, 0))) %>%
    mutate(Pro.T2_01 = 
             ifelse(Pro.T2_total >= 3, 1, 0)) %>%
    mutate(Pro.DWI_01 = 
             ifelse(Pro.DWI_total >= 3, 1, 0)) 
}


#### process data ####

# TB_sup_vs_pro = convert_dates(TB_sup_vs_pro, "ScanDate")

TB_sup_vs_pro = processData_sum_individual_seq(TB_sup_vs_pro)

TB_sup_vs_pro = processData_PIQUALv.2_score(TB_sup_vs_pro)

TB_sup_vs_pro = processData_PIQUALv.2_individual_seq_to_01(TB_sup_vs_pro)


TB_rectal_air = convert_dates(TB_rectal_air, "ScanDate") 


#### pivot longer 

select_metrics = c("No", "AnonID", 
                   "Sup.T2_total", "Sup.DWI_total", "mp.Sup.Overall_3.check", 
                   "Pro.T2_total", "Pro.DWI_total", "mp.Pro.Overall_3.check",
                   "mp.Sup.T2.Pro.DWI.Overall_3.check",
                   "bp.Sup.Overall_3.check", "bp.Pro.Overall_3.check",
                   "bp.Sup.T2.Pro.DWI.Overall_3.check")

TB_sup_vs_pro_long = TB_sup_vs_pro %>%
  select(select_metrics) %>%
  pivot_longer(!(c(No, AnonID)), names_to = "group", values_to = "quality_value")

### join TB_sup_vs_pro with TB_rectal_air

TB_sup_vs_pro = TB_sup_vs_pro %>%
  left_join(TB_rectal_air,
            join_by(No)
  )

### join TB_sup_vs_pro with supine_air_quan


TB_sup_vs_pro = TB_sup_vs_pro %>%
  left_join(supine_air_quan,
            join_by(No)
  )



# change scan date to the correct format
TB_sup_vs_pro <- TB_sup_vs_pro %>%
  mutate(ScanDate.y = case_when(
    !is.na(as.numeric(ScanDate.y)) ~ as.Date(as.numeric(ScanDate.y), origin = "1899-12-30"),
    grepl("^\\d{2}/\\d{2}/\\d{4}$", ScanDate.y) ~ as.Date(ScanDate.y, format = "%d/%m/%Y"),
    TRUE ~ as.Date(NA)  # Assign NA for any values that do not match the expected patterns
  ))

TB_sup_vs_pro %>%
  select(ScanDate.y)

#### Sanity Check ####
## Check overall score (3-point), readers overall score == automatized score

count(TB_sup_vs_pro[TB_sup_vs_pro$mp.Sup.Overall_3 != TB_sup_vs_pro$mp.Sup.Overall_3.check, ])
# TB_sup_vs_pro[TB_sup_vs_pro$Sup.Overall_3 != TB_sup_vs_pro$Sup.Overall_3.check, ]

count(TB_sup_vs_pro[TB_sup_vs_pro$mp.Pro.Overall_3 != TB_sup_vs_pro$mp.Pro.Overall_3.check, ])
# TB_sup_vs_pro[TB_sup_vs_pro$Pro.Overall_3 != TB_sup_vs_pro$Pro.Overall_3.check, ]


TB_sup_vs_pro %>% 
  select(No, AnonID, mp.Sup.Overall_3, mp.Sup.Overall_3.check, mp.Pro.Overall_3, mp.Pro.Overall_3.check) %>%
  print(n = Inf)

summary(TB_sup_vs_pro)


TB_sup_vs_pro %>%
  select(ScanDate.x, ScanDate.y) %>%
  summarise(result = all.equal(ScanDate.x, ScanDate.y))

TB_sup_vs_pro = TB_sup_vs_pro %>%
  select(-ScanDate.y)

TB_sup_vs_pro %>%
  select(Supine.Rectal.Score, SupineAirScore) %>%
  summarise(result = all.equal(Supine.Rectal.Score, SupineAirScore))

TB_sup_vs_pro = TB_sup_vs_pro %>%
  select(-SupineAirScore)


#### Shapiro-Wilk test check normality ####


TB_sup_vs_pro %>% 
  shapiro_test(Sup.T2_total, Sup.DWI_total, mp.Sup.Overall_3.check,
               mp.Sup.T2.Pro.DWI.Overall_3.check,
               Pro.T2_total, Pro.DWI_total, mp.Pro.Overall_3.check,
               bp.Sup.Overall_3.check, bp.Pro.Overall_3.check,
               bp.Sup.T2.Pro.DWI.Overall_3.check) %>%
  mutate(variable = factor(variable, levels = c("Sup.T2_total", "Sup.DWI_total", "mp.Sup.Overall_3.check", 
                                                "Pro.T2_total", "Pro.DWI_total", "mp.Pro.Overall_3.check",
                                                "mp.Sup.T2.Pro.DWI.Overall_3.check",
                                                "bp.Sup.Overall_3.check", "bp.Pro.Overall_3.check",
                                                "bp.Sup.T2.Pro.DWI.Overall_3.check"))) %>%
  arrange(variable)

# all < 0.05 --> non-parametric



#### using the functions ####

Wilcox_test_and_graph = function(data, compare_group){
  
  
  # find min and max possible quality score
  
  max_quality_score = 
    data %>% 
    filter(group %in% compare_group) %>%
    summarise (max(quality_value)) %>%
    pull()
  
  min_quality_score = 
    data %>% 
    filter(group %in% compare_group) %>%
    summarise (min (quality_value)) %>%
    pull()  
  
  
  y.position.self.defined = max_quality_score + 0.5
  
  # Wilcox test
  c.w_TB_compare_group_selected = 
    data %>%
    filter(group %in% compare_group) %>%
    compare_means(quality_value ~ group,
                  method = "wilcox.test", 
                  paired = TRUE,
                  # p.adjust.method = "bonferroni", 
                  .) %>% 
    mutate(y.position = y.position.self.defined) %>%         ###
    mutate(p = round(p, 3)) %>%
    mutate( 
      p = (
        if (p == 0) {p = '<0.001'} 
        else {p = p}
      )) %>%
    mutate(p_display = paste(p, p.signif, sep = " ")) 
  
  
  c.w_TB_compare_group_selected_fig =
    ggplot(data %>%
             filter(group %in% compare_group),
           aes(x = group, y = quality_value)) +
    geom_quasirandom(width=0.1, size = 0.7) +
    scale_x_discrete(limits= compare_group) +
    scale_y_continuous(breaks = seq( min_quality_score, max_quality_score, by=1), 
                       limits=c( min_quality_score-0.5, max_quality_score + 0.7)) +
    theme_classic() + # use black and white theme
    labs(x = "Sequence",                              ###
         y = "Quality score") +                       ###
    stat_pvalue_manual(
      c.w_TB_compare_group_selected , 
      label = "p_display" 
    ) 
  
  print(c.w_TB_compare_group_selected)
  
  print(c.w_TB_compare_group_selected_fig)
  return(c.w_TB_compare_group_selected_fig)
}

# Function to calculate median, Q1, and Q3 for each group and print a summary line
group_summary_stats <- function(data, compare_group) {
  # Calculate median, Q1, and Q3 for each group
  summary <- data %>%
    filter(group %in% compare_group) %>%
    group_by(group) %>%
    summarise(
      median = median(quality_value, na.rm = TRUE),
      Q1 = quantile(quality_value, 0.25, na.rm = TRUE),
      Q3 = quantile(quality_value, 0.75, na.rm = TRUE),
      mean = mean(quality_value, na.rm = TRUE)
    )
  
  # Print the individual summary table without Q3 and median_Q1_Q3 columns
  print(summary)
  
  # Generate and print the combined summary line using median, Q1, and Q3 for each group
  q1_q3_values <- data %>%
    filter(group %in% compare_group) %>%
    group_by(group) %>%
    summarise(
      Q1 = quantile(quality_value, 0.25, na.rm = TRUE),
      Q3 = quantile(quality_value, 0.75, na.rm = TRUE),
      median = median(quality_value, na.rm = TRUE)
    ) %>%
    mutate(median_Q1_Q3 = paste0(median, " [", Q1, ", ", Q3, "]")) %>%
    pull(median_Q1_Q3)
  
  combined_summary <- paste(q1_q3_values, collapse = " vs. ")
  cat("\nCombined summary:", combined_summary, "\n")
}

#### run each comparison by using the function (w/o air threshold) ####
#### "Sup.T2_total", "Pro.T2_total" ####
selected_seq = c("Sup.T2_total", "Pro.T2_total")

c.w_TB_T2_sup_vs_pro_fig = Wilcox_test_and_graph(TB_sup_vs_pro_long, selected_seq)

group_summary_stats(TB_sup_vs_pro_long, selected_seq)

#### save c.w_TB_T2_sup_vs_pro_fig to 300 dpi
tiff("c.w_TB_T2_sup_vs_pro_fig.tiff", units="in", width=5, height=5, res=300)
plot(c.w_TB_T2_sup_vs_pro_fig)        
dev.off()


#### "Sup.DWI_total", "Pro.DWI_total" ####

selected_seq = c("Sup.DWI_total", "Pro.DWI_total")

c.w_TB_DWI_sup_vs_pro_fig = Wilcox_test_and_graph(TB_sup_vs_pro_long, selected_seq)
group_summary_stats(TB_sup_vs_pro_long, selected_seq)

#### save c.w_TB_DWI_sup_vs_pro_fig to 300 dpi
tiff("c.w_TB_DWI_sup_vs_pro_fig.tiff", units="in", width=5, height=5, res=300)
plot(c.w_TB_DWI_sup_vs_pro_fig)        
dev.off()



selected_seq = c("mp.Sup.Overall_3.check", "mp.Pro.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long, selected_seq)

group_summary_stats(TB_sup_vs_pro_long, selected_seq)



selected_seq = c("mp.Sup.Overall_3.check", "mp.Sup.T2.Pro.DWI.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long, selected_seq)

group_summary_stats(TB_sup_vs_pro_long, selected_seq)



selected_seq = c("bp.Sup.Overall_3.check", "bp.Pro.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long, selected_seq)

group_summary_stats(TB_sup_vs_pro_long, selected_seq)



selected_seq = c("bp.Sup.Overall_3.check", "bp.Sup.T2.Pro.DWI.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long, selected_seq)

group_summary_stats(TB_sup_vs_pro_long, selected_seq)



#### air volume < 4.0 cc as a cut off ####

### create tables
Air_volume_less_than_4.0_AnonID = TB_sup_vs_pro %>%
  filter (SupineAirVolume < 4.0) %>%
  select (AnonID, SupineAirVolume)


TB_sup_vs_pro_long_air_volume_less_than_4.0 = TB_sup_vs_pro_long %>%
  filter(AnonID %in% Air_volume_less_than_4.0_AnonID$AnonID)

view(TB_sup_vs_pro_long_air_volume_less_than_4.0)


#### run each comparison by using the function (air volume < 4.0 cc) ####

selected_seq = c("Sup.T2_total", "Pro.T2_total")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)



selected_seq = c("Sup.DWI_total", "Pro.DWI_total")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)



selected_seq = c("mp.Sup.Overall_3.check", "mp.Pro.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)



selected_seq = c("mp.Sup.Overall_3.check", "mp.Sup.T2.Pro.DWI.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)



selected_seq = c("bp.Sup.Overall_3.check", "bp.Pro.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)



selected_seq = c("bp.Sup.Overall_3.check", "bp.Sup.T2.Pro.DWI.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_less_than_4.0, selected_seq)



#### air volume > 4.0 cc as a cut off ####

### create tables
Air_volume_more_than_4.0_AnonID = TB_sup_vs_pro %>%
  filter (SupineAirVolume > 4.0) %>%
  select (AnonID, SupineAirVolume)


TB_sup_vs_pro_long_air_volume_more_than_4.0 = TB_sup_vs_pro_long %>%
  filter(AnonID %in% Air_volume_more_than_4.0_AnonID$AnonID)

view(TB_sup_vs_pro_long_air_volume_more_than_4.0)


#### run each comparison by using the function (air volume > 4.0 cc) ####

selected_seq = c("Sup.T2_total", "Pro.T2_total")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)



selected_seq = c("Sup.DWI_total", "Pro.DWI_total")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)



selected_seq = c("mp.Sup.Overall_3.check", "mp.Pro.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)



selected_seq = c("mp.Sup.Overall_3.check", "mp.Sup.T2.Pro.DWI.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)



selected_seq = c("bp.Sup.Overall_3.check", "bp.Pro.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)



selected_seq = c("bp.Sup.Overall_3.check", "bp.Sup.T2.Pro.DWI.Overall_3.check")

Wilcox_test_and_graph(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_air_volume_more_than_4.0, selected_seq)


##### Buscopan #####
Buscopan_not_given = Buscopan_data %>%
  filter(Buscopan == 0)

compare_group = c("Sup.T2_total", "Pro.T2_total")

TB_sup_vs_pro_long_Buscopan_not_given = 
  TB_sup_vs_pro_long %>%
  filter(AnonID %in% Buscopan_not_given$AnonID)


selected_seq = c("Sup.T2_total", "Pro.T2_total")

c.w_TB_T2_sup_vs_pro_BuscopanNotGiven_fig = 
  Wilcox_test_and_graph(TB_sup_vs_pro_long_Buscopan_not_given, selected_seq)

group_summary_stats(TB_sup_vs_pro_long_Buscopan_not_given, selected_seq)


####### part 3 #######

###### handy functions ######

#### results of individual sequences ####

target_sequence <- function(sequence) {
  # Step 1: Keep original TB and IC datasets, select the user-defined sequence column
  TB_sup_vs_pro_target <- TB_sup_vs_pro %>%
    select("No", "AnonID", "ScanDate", all_of(sequence))
  
  IC_sup_vs_pro_target <- IC_sup_vs_pro %>%
    select("No", "AnonID", "ScanDate", all_of(sequence))
  
  # Step 2: Join the datasets
  TB_IC_sup_vs_pro_combined_target <- TB_sup_vs_pro_target %>%
    inner_join(IC_sup_vs_pro_target, by = c("No", "AnonID", "ScanDate"), suffix = c(".TB", ".IC"))
  
  print(TB_IC_sup_vs_pro_combined_target)  # Display joined data
  
  # Step 3: Generate the contingency table dynamically
  col_tb <- paste0(sequence, ".TB")  # Construct TB column name
  col_ic <- paste0(sequence, ".IC")  # Construct IC column name
  
  # Check if columns exist before proceeding
  if (!(col_tb %in% names(TB_IC_sup_vs_pro_combined_target)) | !(col_ic %in% names(TB_IC_sup_vs_pro_combined_target))) {
    stop("One or both specified columns do not exist in the dataset.")
  }
  
  cat("\n### Contingency Table ###\n")
  print(with(TB_IC_sup_vs_pro_combined_target, table(get(col_tb), get(col_ic))))  # Print contingency table
  
  # Step 4: Remove No, AnonID, and ScanDate for irrCAC calculations
  TB_IC_sup_vs_pro_combined_target_processed <- TB_IC_sup_vs_pro_combined_target %>%
    select(-No, -AnonID, -ScanDate)
  
  print(TB_IC_sup_vs_pro_combined_target_processed)  # Display processed data
  
  # Step 5: Calculate Gwet's AC1 (unweighted)
  cat("\n### Gwet's AC1 (Unweighted) ###\n")
  ac1_unweighted <- gwet.ac1.raw(TB_IC_sup_vs_pro_combined_target_processed, 
                                 weights = "unweighted", categ.labels = NULL, 
                                 conflev = 0.95, N = Inf)
  print(ac1_unweighted$est)
  
  # Step 6: Calculate Gwet's AC1 (Quadratic)
  cat("\n### Gwet's AC1 (Quadratic) ###\n")
  ac1_quadratic <- gwet.ac1.raw(TB_IC_sup_vs_pro_combined_target_processed, 
                                weights = "quadratic", categ.labels = NULL, 
                                conflev = 0.95, N = Inf)
  print(ac1_quadratic$est)
  
  # Return results as a list
  return(list(
    TB_IC_sup_vs_pro_combined_target = TB_IC_sup_vs_pro_combined_target, 
    TB_IC_sup_vs_pro_combined_target_processed = TB_IC_sup_vs_pro_combined_target_processed, 
    contingency_table = with(TB_IC_sup_vs_pro_combined_target, table(get(col_tb), get(col_ic))),
    ac1_unweighted = ac1_unweighted$est, 
    ac1_quadratic = ac1_quadratic$est
  ))
}



#### Use the function:
result_Sup.T2_total <- target_sequence("Sup.T2_total")
result_Sup.DWI_total <- target_sequence("Sup.DWI_total")
result_Pro.T2_total <- target_sequence("Pro.T2_total")
result_Pro.DWI_total <- target_sequence("Pro.DWI_total")

# result_Sup.T2_total$ac1_unweighted
result_Sup.T2_total$ac1_quadratic

# result_Pro.T2_total$ac1_unweighted
result_Pro.T2_total$ac1_quadratic


# result_Sup.DWI_total$ac1_unweighted
result_Sup.DWI_total$ac1_quadratic

# result_Pro.DWI_total$ac1_unweighted
result_Pro.DWI_total $ac1_quadratic
