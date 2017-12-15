Sys.setenv("plotly_username" = "amiatroll")
Sys.setenv("plotly_api_key" = "Kj3wxnvhUFEu1eyO6SgO")
Sys.setenv('MAPBOX_TOKEN' = 'pk.eyJ1IjoiYW1pYXRyb2xsIiwiYSI6ImNqOXQ5MW8yMzNpbzEyd2xnZ281NHpyY2gifQ.oeTfbg1MPdDckUy706ZuIQ')

library(plotly)


# df = read.csv('https://raw.githubusercontent.com/bcdunbar/datasets/master/meteorites_subset.csv')
# 
# p <- df %>%
#   plot_mapbox(lat = ~reclat, lon = ~reclong,
#               split = ~class, size=2,
#               mode = 'scattermapbox', hoverinfo='name') %>%
#   layout(title = 'Meteorites by Class',
#          font = list(color='white'),
#          plot_bgcolor = '#191A1A', paper_bgcolor = '#191A1A',
#          mapbox = list(style = 'dark'),
#          legend = list(orientation = 'h',
#                        font = list(size = 8)),
#          margin = list(l = 25, r = 25,
#                        b = 25, t = 25,
#                        pad = 2))
# 
# # Create a shareable link to your chart
# # Set up API credentials: https://plot.ly/r/getting-started
# chart_link = plotly_POST(p, filename="multiple")
# chart_link

library(dplyr)

# airport locations
air <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/2011_february_us_airport_traffic.csv')

# flights between airports
flights <- read.csv('https://raw.githubusercontent.com/plotly/datasets/master/2011_february_aa_flight_paths.csv')
flights$id <- seq_len(nrow(flights))

p <- plot_mapbox(mode = 'scattermapbox') %>%
  add_markers(
    data = air, x = ~long, y = ~lat, text=~airport, color=I("red"),
    size = ~cnt, hoverinfo = "text", alpha = 0.5) %>%
  add_segments(
    data = group_by(flights, id),
    x = ~start_lon, xend = ~end_lon,
    y = ~start_lat, yend = ~end_lat,
    alpha = 0.3, size = I(1), hoverinfo = "none",
    color=I("red")) %>%
  layout(
    plot_bgcolor = '#191A1A', paper_bgcolor = '#191A1A',
    mapbox = list(style = 'dark',
                  zoom = 1.5,
                  center = list(lat = median(air$lat),
                                lon = median(air$long))),
    margin = list(l = 0, r = 0,
                  b = 0, t = 0,
                  pad = 0),
    showlegend=FALSE)

# Create a shareable link to your chart
# Set up API credentials: https://plot.ly/r/getting-started
chart_link = plotly_POST(p, filename="lines")
chart_link