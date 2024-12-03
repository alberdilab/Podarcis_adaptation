#Hot and dry population
family_summary %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
  group_by(family, Population) %>%
  filter(Population=="Hot_dry")%>%
  filter(time_point=="0_Wild") %>%
  summarise(total_mean=mean(relabun*100, na.rm=T),
            total_sd=sd(relabun*100, na.rm=T))  %>%
  mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
  arrange(-total_mean) %>% 
  dplyr::select(family,total, Population) %>% 
  tt()

family_summary %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
  group_by(family, Population) %>%
  filter(Population=="Cold_wet")%>%
  filter(time_point=="0_Wild") %>%
  summarise(total_mean=mean(relabun*100, na.rm=T),
            total_sd=sd(relabun*100, na.rm=T))  %>%
  mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
  arrange(-total_mean) %>% 
  dplyr::select(family,total, Population) %>% 
  tt()

genus_summary %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
  group_by(genus, Population) %>%
  filter(Population=="Cold_wet")%>%
  filter(time_point=="0_Wild") %>%
  summarise(total_mean=mean(relabun*100, na.rm=T),
            total_sd=sd(relabun*100, na.rm=T))  %>%
  mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
  arrange(-total_mean) %>% 
  dplyr::select(genus,total, Population) %>% 
  tt()

genus_summary %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
  group_by(genus, Population) %>%
  filter(Population=="Hot_dry")%>%
  filter(time_point=="0_Wild") %>%
  summarise(total_mean=mean(relabun*100, na.rm=T),
            total_sd=sd(relabun*100, na.rm=T))  %>%
  mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
  arrange(-total_mean) %>% 
  dplyr::select(genus,total, Population) %>% 
  tt()

genus_summary %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
  group_by(genus, type) %>%
  filter(Population=="Cold_wet")%>%
  filter(time_point=="6_Post-FMT2") %>%
  summarise(total_mean=mean(relabun*100, na.rm=T),
            total_sd=sd(relabun*100, na.rm=T))  %>%
  mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
  arrange(-total_mean) %>% 
  dplyr::select(genus,total, type) %>% 
  tt()


#captivity comparison
phylum_summary %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
  group_by(phylum, Population) %>%
  filter(Population=="Hot_dry")%>%
  filter(time_point=="1_Acclimation") %>%
  summarise(total_mean=mean(relabun*100, na.rm=T),
            total_sd=sd(relabun*100, na.rm=T))  %>%
  mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
  arrange(-total_mean) %>% 
  dplyr::select(phylum,total, Population) %>% 
  tt()

phylum_summary %>%
  left_join(sample_metadata, by = join_by(sample == Tube_code))  %>%
  group_by(phylum, Population) %>%
  filter(Population=="Cold_wet")%>%
  filter(time_point=="1_Acclimation") %>%
  summarise(total_mean=mean(relabun*100, na.rm=T),
            total_sd=sd(relabun*100, na.rm=T))  %>%
  mutate(total=str_c(round(total_mean,2),"±",round(total_sd,2))) %>% 
  arrange(-total_mean) %>% 
  dplyr::select(phylum,total, Population) %>% 
  tt()

wild.counts_plot<-wild.counts %>%
tibble::rownames_to_column(var = "genome")

genome_counts_rel <- wild.counts_plot %>%
  mutate_at(vars(-genome),~./sum(.)) %>%
  column_to_rownames(., "genome")
genome_counts_rel_pa=1*(genome_counts_rel>0)

table_upset_analysis_cont=t(aggregate(t(genome_counts_rel_pa),by=list(sample_metadata_wild$Population),FUN=sum)[,-1])
colnames(table_upset_analysis_cont)=levels(as.factor(sample_metadata_wild$Population))
table_upset_analysis=(table_upset_analysis_cont>0)*1
table_upset_analysis=data.frame(table_upset_analysis)
table_upset_analysis=apply(table_upset_analysis,2,as.integer)
rownames(table_upset_analysis) <- rownames(genome_counts_rel_pa)

locationcolors=c('#008080',"#d57d2c")
upset(as.data.frame(table_upset_analysis),
      keep.order = T,
      sets = rev(c("Cold_wet","Hot_dry")),
      sets.bar.color= rev(locationcolors),
      mb.ratio = c(0.55, 0.45), order.by = "freq")


pha_rep2 %>%
  ggplot(aes(x = time_point, y = normalized, color=time_point, fill=time_point, alpha=0.2 )) +
  geom_boxplot()+
  geom_jitter() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  scale_color_manual(name="Time_point",
                     breaks=c("0","1"),
                     labels=c("Pre","24h"),
                     values=c('#BFA366', "#dec14b")) +
  scale_fill_manual(name="Time_point",
                    breaks=c("0","1"),
                    labels=c("Pre","24h"),
                    values=c('#BFA36650', "#dec14b50"))+
  facet_wrap(~ factor(Type))+
  stat_compare_means(size=3)+
  theme(legend.position="none")+
  labs(x = "Time_point", y= "Normalized measurement (mm/g)")



# pha test correlation

correlation <- cor(pha_test$Mean_value_pre, pha_test$Mean_value_post)
print(correlation) #0.58

##plotting the residuals of the body condition

log_mass <- log(pha_test$Weight) 
log_svl <- log(pha_test$SVL)

model <- lm(log_mass ~ log_svl)
residuals <- resid(model)

residuals_data <- data.frame(log_svl = log_svl, residuals = residuals, Weight=pha_test$Weight)

pha_test_full<-residuals_data %>%
  left_join(pha_test, by="Weight")

pha_test_full %>%
  ggplot(aes(x = residuals, y = difference)) +
  geom_point() +
  geom_smooth(method = lm, formula = y ~ x) +
  facet_wrap(~ Type) +
  stat_poly_eq () +
  labs(x = "Body condition (resd)", y="PHA")

##modelling the residuals with the type of the experiment
model1<-lm(pha_test$difference  ~ pha_test$SVL)
residuals1 <- resid(model1)

residuals_data1 <- data.frame(log_svl = log_svl, residuals = residuals1, Weight=pha_test$Weight)
pha_test_full<-residuals_data1 %>%
  left_join(pha_test_full, by="Weight")

pha_test_full_cold<-pha_test_full  %>%
  filter(Type!="WC")

anova_result <- aov(residuals.x ~ Type, data = pha_test_full_cold)
summary(anova_result)

## Wild domain function

#Biosynthesis
p1 <-merge_gift_wild %>%
  ggplot(aes(x=Population,y=Biosynthesis,color=Population,fill=Population))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  scale_color_manual(name="Population",
                     breaks=c("Cold_wet","Hot_dry"),
                     labels=c("Cold","Hot"),
                     values=c('#008080', "#d57d2c")) +
  scale_fill_manual(name="Population",
                    breaks=c("Cold_wet","Hot_dry"),
                    labels=c("Cold","Hot"),
                    values=c('#00808050', "#d57d2c50")) +
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Population")

#Degradation
p2 <-merge_gift_wild %>%
  ggplot(aes(x=Population,y=Degradation,color=Population,fill=Population))+
  geom_jitter(width = 0.2, size = 1.45, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  scale_color_manual(name="Population",
                     breaks=c("Cold_wet","Hot_dry"),
                     labels=c("Cold","Hot"),
                     values=c('#008080', "#d57d2c")) +
  scale_fill_manual(name="Population",
                    breaks=c("Cold_wet","Hot_dry"),
                    labels=c("Cold","Hot"),
                    values=c('#00808050', "#d57d2c50")) +
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Population")

#Structure
p3 <-merge_gift_wild %>%
  ggplot(aes(x=Population,y=Structure,color=Population,fill=Population))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  scale_color_manual(name="Population",
                     breaks=c("Cold_wet","Hot_dry"),
                     labels=c("Cold","Hot"),
                     values=c('#008080', "#d57d2c")) +
  scale_fill_manual(name="Population",
                    breaks=c("Cold_wet","Hot_dry"),
                    labels=c("Cold","Hot"),
                    values=c('#00808050', "#d57d2c50")) +
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 3, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Population")

#Overall
p4 <-merge_gift_wild %>%
  ggplot(aes(x=Population,y=Overall,color=Population,fill=Population))+
  geom_jitter(width = 0.2, size = 1.5, show.legend = FALSE)+ 
  geom_boxplot(alpha=0.2,outlier.shape = NA, width = 0.5, show.legend = FALSE, coef=0)+
  scale_color_manual(name="Population",
                     breaks=c("Cold_wet","Hot_dry"),
                     labels=c("Cold","Hot"),
                     values=c('#008080', "#d57d2c")) +
  scale_fill_manual(name="Population",
                    breaks=c("Cold_wet","Hot_dry"),
                    labels=c("Cold","Hot"),
                    values=c('#00808050', "#d57d2c50")) +
  stat_compare_means() +
  theme(axis.text.x = element_text(vjust = 0.5, size=10),
        axis.text.y = element_text(size=10),
        axis.title=element_text(size=12,face="bold"),
        axis.text = element_text(face="bold", size=18),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(size = 0.5, linetype = "solid", colour = "black"),
        legend.text = element_text(size=10),
        legend.title = element_text(size=12),
        legend.position="none",
        legend.key.size = unit(1, 'cm'),
        strip.text.x = element_text(size = 12, color = "black", face = "bold"))+
  labs( x = "Population")
```

grid.arrange(arrangeGrob(p1, p2,p3, p4, ncol = 2))
