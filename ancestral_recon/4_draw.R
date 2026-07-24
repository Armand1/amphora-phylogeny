suppressMessages({library(treeio);library(ggtree);library(ape);library(ggplot2);library(dplyr)})
setwd("/Users/armandleroi/Documents/Amphora_phylogeny/phylogenies/2026_07_three_clade_tree")
FIGW<-11; FIGH<-13
trb<-read.beast("runs/srvf38pc/output/simple_run_srvf38PC_mcc.tree"); phy<-as.phylo(trb)
sel<-read.csv("selected_taxa.csv")
meta<-sel|>transmute(label=type,clade=dplyr::recode(sub("_(greek|ital)$","",group),loutro="loutrophoros",pan="panathenaic",outgroup="neck (outgroup)",pelike="pelike"),origin=ifelse(place=="greek_italiote","Italiote","Greek"))
prof<-read.csv("/Users/armandleroi/Documents/Amphora_phylogeny/Making amphora types/05_2026_amphora_297_types/2026_amphora_297_types_profiles_with_metadata.csv")
nc<-read.csv("prep/srvf38_node_curves.csv")
p<-ggtree(trb) %<+% meta
xr<-max(p$data$x); ntip<-Ntip(phy)
asp<-(xr/ntip)*(FIGH/FIGW); sy<-0.80; sx<-sy*asp
norm<-function(d) d|>group_by(gid)|>mutate(h=max(y)-min(y),xn=(x-mean(x))/h,yn=(y-min(y))/h)|>ungroup()
tipxy<-p$data|>filter(isTip)|>select(label,tx=x,ty=y)
tp<-prof|>filter(type%in%phy$tip.label)|>transmute(gid=type,point_order,x,y)|>norm()|>
  rename(label=gid)|>left_join(tipxy,by="label")|>left_join(meta,by="label")|>arrange(label,point_order)|>
  mutate(x=tx+xn*sx,y=ty+(yn-0.5)*sy)
nodexy<-p$data|>filter(!isTip)|>select(node,nx=x,ny=y)
an<-nc|>transmute(gid=node,point_order,x,y)|>norm()|>rename(node=gid)|>
  left_join(nodexy,by="node")|>arrange(node,point_order)|>mutate(x=nx+xn*sx,y=ny+(yn-0.5)*sy)
pal<-c(pelike="#E69F00",loutrophoros="#009E73",panathenaic="#56B4E9",`neck (outgroup)`="#999999")
g<-p+
  geom_polygon(data=an,aes(x,y,group=node),fill="grey82",colour="grey45",linewidth=0.2)+
  geom_polygon(data=tp,aes(x,y,group=label,fill=clade),colour="grey30",linewidth=0.2)+
  scale_fill_manual(values=pal,name="Clade (tips)")+
  coord_cartesian(clip="off")+theme_tree2()+
  labs(title="SRVF-38 tree with reconstructed ancestral profiles",
       subtitle="grey = ancestral shape (BM ASR back-projected through SRVF); coloured = observed tip types")
ggsave("runs/srvf38pc/plots/tree_ancestral_profiles.png",g,width=FIGW,height=FIGH,dpi=135)
cat("done\n")
