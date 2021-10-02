library(openxlsx)
library(dplyr)
library(ggplot2)


#### Read data & clean ####
indata <- bind_rows(
  readxl::read_excel('DataforDavev3.xlsx', sheet = 'Sheet2',  range = 'C2:I130') %>%
    slice(c(-1,-2)) %>%
    transmute(subject=`...1`, ct_e = as.numeric(`...2`), culturable_type='Vero Cells', culturable = (!is.na(`...3`)) | (!is.na(`...4`)),
            ct_sge1=as.numeric(ifelse(culturable, `...3`, `...6`)), ct_sge2 = as.numeric(ifelse(culturable, `...4`, `...7`))),
  readxl::read_excel('DataforDavev3.xlsx', sheet = 'Sheet2',  range = 'K2:Q130') %>%
    slice(c(-1,-2)) %>%
    transmute(subject=`...1`, ct_e = as.numeric(`...2`), culturable_type='TMPRSS2 Cells', culturable = (!is.na(`...3`)) | (!is.na(`...4`)),
              ct_sge1=as.numeric(ifelse(culturable, `...3`, `...6`)), ct_sge2 = as.numeric(ifelse(culturable, `...4`, `...7`)))
) 


extra_vals <- readxl::read_excel('DataforDavev3.xlsx', sheet = 'Sheet2',  range = 'S1:W106') %>%
  transmute(subject=`Code`, ct_e = `E (orig)`,
            neg_e = as.numeric(case_when(trimws(`Neg E`) == 'NDET' ~ '45', 
                                         trimws(`Neg E`) == '--' ~ NA_character_, 
                                         TRUE ~ `Neg E`)), 
            ct_e_repeat = as.numeric(case_when(trimws(`E (new)`) == 'NDET' ~ '45', 
                                               trimws(`E (new)`) == '--' ~ NA_character_, 
                                               TRUE ~ `E (new)`)))


indata_formatted <- indata %>% 
  mutate(subject=sub('\\*|\\^', '', subject)) %>%
  mutate_at(vars(ct_e, ct_sge1, ct_sge2), function(x) case_when(x>=40 ~ 45, TRUE ~ x)) %>%
  reshape2::dcast(subject + ct_e + ct_sge1 + ct_sge2 ~ paste('Culturable in', culturable_type), value.var = 'culturable') %>%
  as_tibble() %>%
  left_join(extra_vals %>% select(-ct_e), by = c("subject"))

#### Build ROC curve data ####
graph_data <- bind_rows(
  indata_formatted %>%
    mutate(culturable=`Culturable in TMPRSS2 Cells`) %>%
    arrange(ct_e, culturable) %>%
    filter(!is.na(ct_e)) %>%
    transmute(culturables = cumsum(culturable==TRUE),
              tot_culturables = sum(culturable==TRUE),
              nonculturables = cumsum(culturable==FALSE),
              tot_nonculturables = sum(culturable==FALSE),
              true_pos=cumsum(culturable==TRUE)/sum(culturable==TRUE),
              false_pos=cumsum(culturable==FALSE)/sum(culturable==FALSE),
              thresh=ct_e,
              method='ct_e',
              cell_type='TMPRSS2'),
  indata_formatted %>%
    mutate(culturable=`Culturable in TMPRSS2 Cells`) %>%
    arrange(neg_e, culturable) %>%
    filter(!is.na(neg_e)) %>%
    transmute(culturables = cumsum(culturable==TRUE),
              tot_culturables = sum(culturable==TRUE),
              nonculturables = cumsum(culturable==FALSE),
              tot_nonculturables = sum(culturable==FALSE),
              true_pos=cumsum(culturable==TRUE)/sum(culturable==TRUE),
              false_pos=cumsum(culturable==FALSE)/sum(culturable==FALSE),
              thresh=neg_e,
              method='neg_e',
              cell_type='TMPRSS2'),
  indata_formatted %>%
    mutate(culturable=`Culturable in TMPRSS2 Cells`) %>%
    arrange(ct_sge1, culturable) %>%
    filter(!is.na(ct_sge1)) %>%
    transmute(culturables = cumsum(culturable==TRUE),
              tot_culturables = sum(culturable==TRUE),
              nonculturables = cumsum(culturable==FALSE),
              tot_nonculturables = sum(culturable==FALSE),
              true_pos=cumsum(culturable==TRUE)/sum(culturable==TRUE),
              false_pos=cumsum(culturable==FALSE)/sum(culturable==FALSE),
              thresh=ct_sge1,
              method='ct_sge1',
              cell_type='TMPRSS2'),
  indata_formatted %>%
    mutate(culturable=`Culturable in TMPRSS2 Cells`) %>%
    arrange(ct_sge2, culturable) %>%
    filter(!is.na(ct_sge2)) %>%
    transmute(culturables = cumsum(culturable==TRUE),
              tot_culturables = sum(culturable==TRUE),
              nonculturables = cumsum(culturable==FALSE),
              tot_nonculturables = sum(culturable==FALSE),
              true_pos=cumsum(culturable==TRUE)/sum(culturable==TRUE),
              false_pos=cumsum(culturable==FALSE)/sum(culturable==FALSE),
              thresh=ct_sge2,
              method='ct_sge2',
              cell_type='TMPRSS2')
)



#### Save ROC data ####
wb <- openxlsx::createWorkbook(creator = 'Dave Shirley')
openxlsx::addWorksheet(wb, 'ROC Data')
openxlsx::writeDataTable(wb, 'ROC Data', graph_data,
                         tableStyle = 'none', 
                         headerStyle = openxlsx::createStyle(textDecoration = 'bold', border='bottom', borderStyle = 'medium'))
openxlsx::saveWorkbook(wb, 'ROC Data.xlsx', overwrite = TRUE)


#### make ROC graph ####
graph_data %>%
  filter(!method %in% c('ct_e_repeat', 'ld_e_sge2')) %>%
  filter(cell_type=='TMPRSS2') %>%
  mutate(method=forcats::fct_recode(method, 
                                    `CT(E)` = 'ct_e', `CT(sgE1)`='ct_sge1', `CT(sgE2)`='ct_sge2',
                                    `CT(negative strand E)`='neg_e')) %>%
  group_by(method, cell_type) %>%
  mutate(label=case_when(
    seq_len(n())==which.max(true_pos>=0.95) ~ thresh, TRUE ~ NA_real_)) %>% #first TRUE
  ungroup() %>%
  ggplot(aes(x=false_pos, y=true_pos,
             color=method, group=interaction(method, thresh<40))) +
  geom_line(aes(linetype=forcats::fct_rev(factor(!thresh>40)))) +
  ggrepel::geom_label_repel(aes(label=round(label, digits = 2))) +
  geom_hline(yintercept = 0.95, linetype=2) +
  theme_bw() + 
  facet_wrap(~cell_type) +
  coord_equal() +
  scale_x_continuous(labels=scales::percent) +
  scale_y_continuous(labels=scales::percent) +
  labs(x='False positive rate: Considered positive based on CT but negative culture
       FP >> quarantine but not actually infectious',
       y='True positive rate: Considered positive based on CT and positive culture
       1-TP=FN >> no quarantine but infectious',
       linetype='CT within LoD',
       color='Marker')
ggsave('graphs/ROC_TMPRSS2_only.pdf', height = 7, width=9)


## truth = culture
## test = ct thresh
## all samples are positive for E, but not necessarily infectious at time of sampling
## true positive = test +ve | culture +ve >> correct recommendation, quarantine, safe
## false positive = test +ve | culture -ve >> wrong recommendation, quarantine, safe (overly conservative)
## true negative = test -ve | culture -ve >> correct recommendation, don't quarantine, safe
## false negative = test -ve | culture +ve >> wrong recommendation, don't quarantine, NOT SAFE

# false negative = 1- true positive