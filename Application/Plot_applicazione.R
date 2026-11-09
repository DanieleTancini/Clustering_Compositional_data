library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)
library(scales)
library(readxl)
library(spdep)


load("K3_apl_sum_exch.RData")
rm(list = setdiff(ls(), c("U_est")))
Dataset_compositional = read_excel("Dati_provinciali.xlsx")

Data = Dataset_compositional[,10:15]
Territorio = Dataset_compositional$TERRITORIO[1:107]
Y = array(NA, dim = c(6,107,4))

for(t in 1:4){
  for(i in 1:107){
    Y[,i,t] = as.numeric(Data[107*(t-1)+i,])
  }
}

dim(Y)
Y[3, , ] = Y[3,,]+Y[6,,]
Y = Y[-6,,]

sectors = c(
  "Agriculture",
  "Industry",
  "Commerce and transport",
  "Financial services",
  "Other services"
)

#convert array to data frame
df = expand.grid(sector = sectors, province = 1:dim(Y)[2], year = 2020:2023)

df$share = as.vector(Y)

df$sector = factor(df$sector, levels = sectors)

ggplot(df, aes(x = factor(province), y = share, fill = sector)) +
  geom_col(width = 0.85) +
  facet_wrap(~year, ncol = 2) +
  scale_y_continuous(labels = percent_format()) +
  labs(
    x = "Province",
    y = "Share of total value added",
    fill = ""
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),  # Remove province numbers
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    strip.text = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 13),
    panel.grid.major.x = element_blank(),
    legend.position = "bottom"
  )

#------------------------

Dataset_compositional = read_excel("Dati_provinciali.xlsx")

Data = Dataset_compositional[,10:15]
Territorio = Dataset_compositional$TERRITORIO[1:107]
Y = array(NA, dim = c(6,107,4))

for(t in 1:4){
  for(i in 1:107){
    Y[,i,t] = as.numeric(Data[107*(t-1)+i,])
  }
}

dim(Y)
Y[3, , ] = Y[3,,]+Y[6,,]
Y = Y[-6,,]

# Read shapefile
prov = st_read("ProvCM01012025_g_WGS84.shp")

# Check data
names(prov)
plot(st_geometry(prov))

nb = poly2nb(prov, queen = TRUE)
W = nb2mat(nb, style = "B", zero.policy = TRUE)
rownames(W) = prov$DEN_UTS
colnames(W) = prov$DEN_UTS

#Uniforma i nomi di W a quelli di Territorio
nm = rownames(W)

nm[nm == "Aosta"] = "Valle d'Aosta"
nm[nm == "Forli'-Cesena"] = "Forlì"
nm[nm == "Massa Carrara"] = "Massa-Carrara"
nm[nm == "Pesaro e Urbino"] = "Pesaro-Urbino"
nm[nm == "Reggio di Calabria"] = "Reggio Calabria"
nm[nm == "Reggio nell'Emilia"] = "Reggio Emilia"
nm[nm == "Verbano-Cusio-Ossola"] = "Verbania"

rownames(W) = nm
colnames(W) = nm

prov$DEN_UTS = nm

#--------

order_ind = order(Territorio)
Territorio[order_ind]
Y_fin = Y[,order_ind,]
Y = Y[,order_ind,]
prov = prov[order(rownames(W)),]
W = W[order(rownames(W)), order(rownames(W))]

#check
Territorio[order_ind] == rownames(W)

prov$DEN_UTS == Territorio[order_ind]

#------------------

library(sf)
library(ggplot2)
library(dplyr)
library(tidyr)

prov_long = bind_rows(lapply(1:4, function(t){
    prov %>%
      mutate(
        State = factor(
          U_est[, t],
          levels = 1:3,
          labels = c("State 1", "State 2", "State 3")
        ),
        Year = factor(2020+t-1)
      )
  })
)

ggplot(prov_long) +
  geom_sf(aes(fill = State), color = "white", linewidth = 0.15) +
  scale_fill_manual(
    values = c(
      "State 1" = "#B2182B",
      "State 2" = "#EF8A62",
      "State 3" = "#FDCC8A"
    ),
    name = "",
    drop = FALSE
  ) +
  facet_wrap(~ Year, ncol = 2) +
  coord_sf(datum = NA) +
  theme_void() +
  theme(
    strip.text = element_text(face = "bold", size = 18),
    legend.position = "bottom",
    legend.title = element_text(face = "bold"),
    legend.text = element_text(size = 18),
    plot.margin = margin(5, 5, 5, 5)
  )
