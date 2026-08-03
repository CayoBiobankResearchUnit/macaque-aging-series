# The theme I use to make consistent font sizes.
require(ggplot2); require(dplyr); require(tidyr); require(cowplot)

# cool https://joeystanley.com/blog/custom-themes-in-ggplot2
pub_theme <- function() {
  theme_cowplot(font_size = 8) %+replace%
    theme(
      strip.background = element_blank(),
      text = element_text(size=8),
      strip.text = element_text(size=7, margin=margin(b=5, t=5)),
      axis.text.x = element_text(size=8,angle = 90),
      axis.text.y = element_text(size=8),
      legend.text = element_text(size=8),
      legend.title = element_text(size=8,hjust=0.5),
      axis.line = element_line(linewidth = 0.4),
      axis.title = element_text(size=8),
      legend.position = 'bottom', 
      plot.tag = element_text(size = 10, face='bold'),
      plot.title = element_text(hjust = 0.5, face='plain', size=9,
                                margin=margin(b=10))
    )
}
theme_set(pub_theme())

#group color
library(dplyr)
color_map = read.csv('/data/CEM/smacklab/team_scratch/wyang114/data/colors.txt', sep='\t',header=FALSE)
colnames(color_map) = c('identity','color')
ct_color = color_map[c(44:69),]
ct_color = ct_color%>%
  mutate(group = case_when(identity == 'GABAergic neurons' ~ 'gaba',
                           identity == 'astrocytes' ~ 'astrocytes',
                           identity == 'basket cells' ~'basket',
                           identity == 'cerebellar neurons' ~'cerebellar',
                           identity == 'ependymal cells' ~ 'ependymal',
                           identity =='glutamatergic neurons' ~ 'glut',
                           identity =='medium spiny neurons' ~ 'msn',
                           identity =='microglia' ~ 'microglia',
                           identity == 'midbrain neurons' ~ 'midbrain',
                           identity == 'oligodendrocytes' ~ 'oligodendrocytes',
                           identity == 'oligodendrocyte precursor cells' ~ 'oligodendrocytes',
                           identity == 'vascular cells' ~ 'vascular'))
region_color = color_map[c(14:43),]
