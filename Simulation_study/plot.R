#plot pdf 12, 12
library(plotly)
load("Simulation_1_Dirichlet_250_edges.RData")


library(igraph)

g = graph_from_adjacency_matrix(Com, mode = "undirected", diag = FALSE)

plot(g, vertex.size = 12, vertex.label = NA, edge.width = 2)

#convert the 4-part composition to 3D coordinates
#y1 = 1 - y2 - y3 - y4
x_ = y[1,]
y_ = y[2,]
z_ = y[3,]

dat = data.frame(x_,y_,z_,c(u_real$U_est[,,10000]))

#pdf 12x12
library(ggtern)

ggtern(data = dat, aes(x = x_, y = y_, z = z_, color = factor(c(u_real$U_est[,,10000])))) +
  geom_point(alpha = 1, size = 4) +
  labs(x = expression(y[1]), y = expression(y[2]), z = expression(y[3]), color = "State") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )


#-----------------------

#convert the 4-part composition to 3D coordinates
x_ = y[1,]
y_ = y[2,]
z_ = y[4,]

dat = data.frame(x_,y_,z_,c(u_real$U_est[,,10000]))

#pdf 12x12

ggtern(data = dat, aes(x = x_, y = y_, z = z_,color = factor(c(u_real$U_est[,,10000])))) +
  geom_point(alpha = 1, size = 4) +
  labs(
    x = expression(y[1]),
    y = expression(y[2]),
    z = expression(y[4]),
    color = "State"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )

#-----------------------

#convert the 4-part composition to 3D coordinates
x_ = y[2,]
y_ = y[3,]
z_ = y[4,]

dat = data.frame(x_,y_,z_,c(u_real$U_est[,,10000]))

#pdf 12x12

ggtern(data = dat, aes(x = x_, y = y_, z = z_, color = factor(c(u_real$U_est[,,10000])))) +
  geom_point(alpha = 1, size = 4) +
  labs(
    x = expression(y[2]),
    y = expression(y[3]),
    z = expression(y[4]),
    color = "State"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )

#------------------------------------------------------------

library(plotly)
load("Simulation_3_Dirichlet_250_edges.RData")

library(igraph)

g = graph_from_adjacency_matrix(Com, mode = "undirected", diag = FALSE)

plot(g, vertex.size = 12, vertex.label = NA, edge.width = 2)

#convert the 4-part composition to 3D coordinates
#y1 = 1 - y2 - y3 - y4
x_ = y[1,]
y_ = y[2,]
z_ = y[3,]

dat = data.frame(x_,y_,z_,c(u_real$U_est[,,10000]))

#pdf 12x12

ggtern(data = dat, aes( x = x_,y = y_, z = z_, color = factor(c(u_real$U_est[,,10000])))) +
  geom_point(alpha = 1, size = 4) +
  labs(
    x = expression(y[1]),
    y = expression(y[2]),
    z = expression(y[3]),
    color = "State"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )

#-----------------------

#convert the 4-part composition to 3D coordinates

x_ = y[1,]
y_ = y[2,]
z_ = y[4,]

dat = data.frame(x_,y_,z_,c(u_real$U_est[,,10000]))

#pdf 12x12

ggtern(data = dat, aes(x = x_, y = y_, z = z_, color = factor(c(u_real$U_est[,,10000])))) +
  geom_point(alpha = 1, size = 4) +
  labs(
    x = expression(y[1]),
    y = expression(y[2]),
    z = expression(y[4]),
    color = "State"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )

#-----------------------

#convert the 4-part composition to 3D coordinates
x_ = y[2,]
y_ = y[3,]
z_ = y[4,]

dat = data.frame(x_,y_,z_,c(u_real$U_est[,,10000]))


#pdf 12x12

ggtern(data = dat, aes(x = x_, y = y_, z = z_, color = factor(c(u_real$U_est[,,10000])))) +
  geom_point(alpha = 1, size = 4) +
  labs(
    x = expression(y[2]),
    y = expression(y[3]),
    z = expression(y[4]),
    color = "State"
  ) +
  theme_bw() +
  theme(
    axis.title = element_text(size = 30),
    axis.text = element_text(size = 30),
    legend.title = element_text(size = 30),
    legend.text = element_text(size = 30)
  )



