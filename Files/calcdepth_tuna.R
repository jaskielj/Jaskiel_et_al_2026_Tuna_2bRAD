#install.packages("tidyverse") ## Install tidyverse if you don't have it installed yet
library(tidyverse)

setwd("/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/")
basedir <- "/projectnb/mullenl/jaskiel/2bRAD_fastq/Combined_2bRAD_analysis/" # Make sure to edit this to match your $BASEDIR
bam_list <- read_lines(paste0(basedir,"depthlist"))

for (i in 1:length(bam_list)){
  bamfile = bam_list[i]
  # Compute depth stats
  depth <- read_tsv(paste0(basedir, bamfile), col_names = F)$X1
  mean_depth <- mean(depth)
  sd_depth <- sd(depth)
  mean_depth_nonzero <- mean(depth[depth > 0])
  mean_depth_within2sd <- mean(depth[depth < mean_depth + 2 * sd_depth])
  median <- median(depth)
  presence <- as.logical(depth)
  proportion_of_reference_covered <- mean(presence)
  output_temp <- tibble(bamfile, mean_depth, sd_depth, mean_depth_nonzero, mean_depth_within2sd, median, proportion_of_reference_covered)
  
  # Bind stats into dataframe and store sample-specific per base depth and presence data
  if (i==1){
    output <- output_temp
    total_depth <- depth
    total_presence <- presence
  } else {
    output <- bind_rows(output, output_temp)
    total_depth <- total_depth + depth
    total_presence <- total_presence + presence
  }
}

calcdepth_output_tuna = output %>%
  mutate(across(where(is.numeric), round, 3))
write.csv(calcdepth_output_tuna, file = "calcdepth_tuna_all.csv")

save.image(file="tuna_depth_all.RData")
load(file="tuna_depth_all.RData")

#plot output
p<-ggplot(data=output, aes(x=bamfile, y=mean_depth)) +
  geom_bar(stat="identity")
p
p<-ggplot(data=output, aes(x=bamfile, y=proportion_of_reference_covered)) +
  geom_bar(stat="identity")
p
p<-ggplot(data=output, aes(x=bamfile, y=mean_depth_nonzero)) +
  geom_bar(stat="identity")
p
mean(output$mean_depth_nonzero)
#17.1854

# Plot the depth distribution (this may take a few minutes to run)
tibble(total_depth = total_depth, position = 1:length(total_depth))  %>%
  ggplot(aes(x = position, y = total_depth)) +
  geom_point(size = 0.1)

#Rerunning without geom_point options because that's what the error seems to be
tibble(total_depth = total_depth, position = 1:length(total_depth))  %>%
  ggplot(aes(x = position, y = total_depth)) +
  geom_point()

# Total depth per site across all individuals 
total_depth_summary <- count(tibble(total_depth = total_depth), total_depth)
total_presence_summary <- count(tibble(total_presence = total_presence), total_presence)
total_depth_summary %>%
  ggplot(aes(x = total_depth, y = n)) +
  geom_point()
total_depth_summary %>%
  ggplot(aes(x = total_depth, y = n)) +
  geom_point() +
  coord_cartesian(xlim=c(NA, 20))
total_presence_summary %>%
  ggplot(aes(x = total_presence, y = n)) +
  geom_col()
