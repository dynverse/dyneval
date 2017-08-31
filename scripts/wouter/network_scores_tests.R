net1 <- dyngen:::generate_toy_milestone_network("linear")
net2 <- dyngen:::generate_toy_milestone_network("cycle")


### NETCOM ICS
# Does not work as expected, gives the same score for linear and cycles
net_adj1 <- net1 %>% igraph::graph_from_data_frame() %>% igraph::as_adj() %>% as.matrix()
net_adj2 <- net2 %>% igraph::graph_from_data_frame() %>% igraph::as_adj() %>% as.matrix()

alignment <- netcom::align(
  net_adj1,
  net_adj2
)

alignment$score

netcom::ics(net_adj1, net_adj2, alignment$alignment)



## NETDIST NETDIS
# Does not work for small networks (gives NaN because no large graphlets)
graphlet_size <- 4
neighbourhood_size <- 10
centred_graphlet_counts1 <- netdist::netdis_centred_graphlet_counts(net1 %>% igraph::graph_from_data_frame(), graphlet_size, neighbourhood_size)
centred_graphlet_counts2 <- netdist::netdis_centred_graphlet_counts(net2 %>% igraph::graph_from_data_frame(), graphlet_size, neighbourhood_size)

netdist::netdis(centred_graphlet_counts1, centred_graphlet_counts2, graphlet_size = graphlet_size)


## NETDIST NET_EMD
gdd1 <- netdist::gdd(net1 %>% igraph::graph_from_data_frame())
gdd2 <- netdist::gdd(net2 %>% igraph::graph_from_data_frame())
netdist::net_emd(gdd1, gdd2)
