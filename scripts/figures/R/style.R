# Original lineage hues retained with increased saturation; clock palette revised.
LINEAGE_COL <- c(Nichols='#3264CF', SS14='#D84A5B', 'Other/unknown'='grey70')
CONTINENT_COL <- c(Africa='#fec44f',Asia='#de2d26',Caribbean='#9f79ee',Europe='#2171b5',
                  'N.America'='#74c476',Oceania='#c51b8a','S.America'='#ec7014','Other/unknown'='grey70')
FOCAL_COL <- c(A506T='#B2182B',A506V='#EF8A62',G708S='#EF8A62',M625V='#EF8A62',P564L='#2166AC',S394R='#2166AC')
CLOCK_COL <- c(Strict='#8B73AE',Relaxed='#607D8B')
STATE_COL <- c(A='#F5F5F5',T='#B2182B',V='#EF8A62',M='#756BB1',Missing='#BBBBBB')
FONT <- Sys.getenv('PBP_FIGURE_FONT','Nimbus Sans')
# The manuscript source defines standard.textsize=9 and text.size.within=2.4.
BASE_SIZE <- 9
common_theme <- function() {
 theme_bw(base_size=BASE_SIZE,base_family=FONT)+theme(
  text=element_text(size=9),axis.text=element_text(size=8,colour='grey20'),
  axis.title=element_text(size=9),panel.border=element_rect(colour='grey35',linewidth=.35),
  panel.grid.major=element_line(colour='grey92',linewidth=.25),panel.grid.minor=element_blank(),
  strip.background=element_rect(fill='grey85',colour='grey35',linewidth=.3),
  strip.text=element_text(size=9,face='bold'),plot.title=element_text(size=10,face='bold',hjust=.5),
  plot.subtitle=element_text(size=8,colour='grey35'),
  legend.text=element_text(size=8),legend.title=element_text(size=8),
  legend.key.size=grid::unit(.55,'line'),legend.position='bottom',
  plot.margin=margin(5,6,5,6),plot.tag=element_text(size=11,face='bold'))
}
