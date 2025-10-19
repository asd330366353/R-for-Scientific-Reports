# 创建与您提供样式完全一致的流程图
rm(list = ls())
gc()
setwd("C:/Users/唐/Desktop/Scientific Reports")


library(DiagrammeR)

library(DiagrammeRsvg)
library(rsvg)
flowchart <- grViz("
digraph exact_replica {
  
  graph [layout = dot,
         rankdir = TB,
         nodesep = 0.3,
         ranksep = 0.5]
  
  node [shape = rectangle,
        style = solid,
        fontname = 'Arial',
        fontsize = 11,
        width = 3.5,
        height = 0.8,
        fixedsize = false,
        color = 'black',
        fontcolor = 'black']
  
  edge [color = 'black',
        arrowhead = normal,
        arrowsize = 0.7]
  
  # 修正后的节点标签
  A [label = 'Participants included at baseline from CHARLS in 2011 (n=17,705)']
  B [label = 'Exclusion: Age<45 (n=391)']
  C [label = '17,314 participants for further screening']
  D [label = 'Exclusion: Participants were already diagnosed with stroke (n=180), hypertension (n=4,450), heart problems (n=1,034), and diabetes (n=411) at baseline.']
  E [label = '11,239 participants for further screening']
  F [label = 'Exclusion: Missing follow-up data on cardiometabolic diseases in 2020 (n=4,768).']
  G [label = '6,471 participants included in the study']
  
  A -> B -> C -> D -> E -> F -> G
}
")
# 显示流程图
flowchart
flowchart %>% 
  export_svg() %>% 
  charToRaw() %>% 
  rsvg_png("flowchart_highres.png", width = 2000, height = 1200)
