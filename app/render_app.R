# SETUP -------------------------------------------------------------------

setwd("C:/Users/bento/gis630")
load(file="artifacts/summary_img.RData")
library(data.table)
library(dplyr)
library(purrr)
library(terra)
library(caret)
library(sf)
library(spatstat)
library(shiny)
library(ggspatial)
source("R/load_preprocessed_data.R")
source("R/load_bioclim_masks.R")
setwd("C:/Users/bento/gis630/app")


get.mtype.name <- function(m) case_when(
  m == "ipp_glm_mpl_2" ~ "ipp",
  m == "maxent" ~ "maxent",
  m == "logreg" ~ "logistic regression",
  m == "knn" ~ "knn",
  m == "tree" ~ "classification tree",
  m == "randomforest" ~ "random forest",
  m == "xgboost" ~ "xgboost",
  T ~ "")

pred.rast.path <- "../artifacts/prediction_rasters"
pred.rasters <- c("ipp_glm_mpl_2", "maxent", 
                  "logreg", "knn", "tree", 
                  "randomforest", "xgboost") %>%
  set_names() %>%
  purrr::map(function(m) {
    m.type <- get.mtype.name(m)
    species %>%
      set_names %>%
      purrr::map(function(spec) {
        states %>%
          set_names %>%
          purrr::map(function(st) {
            r.path <- file.path(pred.rast.path, 
                                paste0(spec, "_", st, "_", 
                                       m.type, ".tif"))
            if (!file.exists(r.path)) {
              pred.raster <- NULL
            } else {
              pred.raster <- terra::rast(r.path)
            }
            pred.raster
          })
      })
  })
names(pred.rasters) <- get.mtype.name(names(pred.rasters))

spec.state.mods <- expand.grid(species=species, 
                               state=states, 
                               model.type=names(pred.rasters),
                               stringsAsFactors=F) %>%
  as.data.table() %>%
  .[!(model.type == "ipp" & species == "Sanderling" & state == "CO")]


# Function to plot a raster with ggplot
plot.raster.gg <- function(r, fill=NULL, scale.c=T, 
                           legend.title="probability",
                           title=NULL) {
  df <- terra::as.data.frame(r, xy = T)
  if (is.null(fill)) {
    fill <- names(r)[[1]]
  } else {
    p <- ggplot(df, aes(x = x, y=y, fill=get(fill)))
  }
   p <- p + 
    geom_raster() +
    coord_cartesian() +
    labs(fill=legend.title)
  if (scale.c) p <- p + scale_fill_viridis_c()
  if (!is.null(title)) p <- p + labs(title=title)
  p
}

plt.width.height <- function(st) {
  width <- case_when(st == "VT"~6.5*2,
                     st == "NC"~16*2,
                     st == "CO"~12*2,
                     st == "OR"~10*2,
                     T~8*2)
  height <- case_when(st == "VT"~7,
                      st == "NC"~7,
                      st == "CO"~9,
                      st == "OR"~8,
                      T~8)
  list(
    w=width,
    h=height
  )
}

get.model.type.formatted <- function(model.type) {
  case_when(model.type == "ipp" ~ "IPP",
            model.type == "maxent" ~ "MaxEnt",
            model.type == "logistic regression" ~ "Logistic Regression",
            model.type == "knn" ~ "KNN",
            model.type == "classification tree" ~ "Classification Tree",
            model.type == "random forest" ~ "Random Forest",
            model.type == "xgboost" ~ "XGBoost",
            T ~ "")
}

plot.to.svg <- function(p, st, encode=F, plt.dir="plots", plt.name="plt", wh=NULL) {
  # Create a temporary file
  tmpfile <- tempfile(fileext = ".svg")
  # Save the ggplot to this file
  if (is.null(wh)) {
    wh <- plt.width.height(st)
  }
  ggsave(filename = tmpfile, p, width=wh$w, height=wh$h)
  if (encode) {
    # Use img() function from htmltools to display it within a div
    img.src <- paste0("data:image/svg+xml;base64,", base64enc::base64encode(tmpfile))
  } else {
    img.src <- file.path(plt.dir, paste0(plt.name, ".svg"))
    file.copy(tmpfile, img.src, overwrite = T)
  }
  file.remove(tmpfile)
  img.src
}

# Get data into normal data frame format (not `spatstat`)
data <- states %>%
  set_names() %>%
  purrr::map(function(st) {
    # Get raster by state
    r <- rasters[[st]]
    species %>% 
      set_names() %>%
      purrr::map(function(spec) {
        # Load `spatstat` quad data
        Q <- readRDS(file.path("../artifacts", "train_spatstat_Q_2",
                               paste0(st, "_", spec, "_Q.rds")))
        Q.test <- readRDS(file.path("../artifacts", "test_spatstat_Q_2",
                                    paste0(st, "_", spec, "_Q.rds")))
        # Load presence/absence data
        pres.train <- data.table(
          x=Q$data$x, 
          y=Q$data$y,
          presence=factor(T, levels=c(F,T), 
                          labels=c("Absence", "Presence")))
        abs.train <- data.table(
          x=Q$dummy$x, 
          y=Q$dummy$y, 
          presence=factor(F, levels=c(F,T),
                          labels=c("Absence", "Presence")))
        pres.test <- data.table(
          x=Q.test$data$x, 
          y=Q.test$data$y,
          presence=factor(T, levels=c(F,T),
                          labels=c("Absence", "Presence")))
        abs.test <- data.table(
          x=Q.test$dummy$x, 
          y=Q.test$dummy$y,
          presence=factor(F, levels=c(F,T),
                          labels=c("Absence", "Presence")))
        purrr::walk(names(r), function(n) {
          pres.train[, (n) := terra::extract(r[[n]], 
                                             cbind(pres.train$x,
                                                   pres.train$y))]
          abs.train[, (n) := terra::extract(r[[n]], 
                                            cbind(abs.train$x,
                                                  abs.train$y))]
          pres.test[, (n) := terra::extract(r[[n]], 
                                            cbind(pres.test$x,
                                                  pres.test$y))]
          abs.test[, (n) := terra::extract(r[[n]], 
                                           cbind(abs.test$x,
                                                 abs.test$y))]
        })
        list(
          train=data.table::rbindlist(l=list(pres.train, abs.train)) %>%
            na.omit(),
          test=data.table::rbindlist(l=list(pres.test, abs.test)) %>%
            na.omit()
        )
      })
  })


data.sf <- purrr::map_df(1:nrow(spec.state), function(i) {
  spec <- spec.state[i,]$species
  st <- spec.state[i,]$state
  d <- rbindlist(list(data[[st]][[spec]]$train[, .(x, y, presence)],
                      data[[st]][[spec]]$test[, .(x, y, presence)]))
  sf::st_as_sf(d, coords = c("x", "y"), 
               crs=4326, sf_column_name = "geometry") %>%
    mutate(species=spec, state=st)
})

states.sf <- states %>%
  set_names() %>%
  purrr::map(function(st) {
    sf::read_sf(paste0("../data/US_State_Boundaries/", st, "_State_Boundaries.shp")) %>%
      transmute(state=STATE_ABBR) %>%
      sf::st_transform(crs = crs(data.sf))
  })

# OBSERVATION PLOTS -------------------------------------------------------

get.pres.abs.plt <- function(data.sf, states.sf, st, spec, 
                             pres="Presence", title=NULL) {
  sdf <- data.sf %>% filter(species == spec & state == st)
  if (pres != "both") {
    sdf <- sdf %>% filter(presence == pres)
  }
  st.sdf <- states.sf[[st]]
  extent <- sf::st_bbox(st.sdf)
  p <- ggplot(data=st.sdf) +
    geom_sf(fill="lightblue") + 
    geom_sf(fill=NA, color="darkgray") +
    coord_sf(xlim = c(extent[[1]] - .05, extent[[3]] + .05), 
             ylim = c(extent[[2]] - .05, extent[[4]] + .05),
             crs=crs(st.sdf))
  if (pres == "both") {
    p <- p + geom_sf(data=sdf, aes(color=presence, shape=presence)) + 
      scale_shape_manual(values=c("Presence" = "+", "Absence" = "o")) + 
      scale_color_manual(values=c("Presence" = "black", "Absence" = "darkred")) +
      labs(color = "Observation Type", shape = "Observation Type")
  } else if (pres == "Presence") {
    p <- p + geom_sf(data=sdf, color="black", shape="+")
  } else if (pres == "Absence") {
    p <- p + geom_sf(data=sdf, color="darkred", shape="o")
  }
  if (!is.null(title)) {
    p <- p + labs(title=title)
  }
  p
}

refresh.obs.plts <- T
obs.list <- purrr::map(1:nrow(spec.state), function(i) {
  st <- spec.state[i,]$state
  spec <- spec.state[i,]$species
  width <- case_when(st == "VT"~7.5,
                     st == "NC"~14,
                     st == "CO"~10.5,
                     st == "OR"~10.6,
                     T~8) / 1.35
  height <- case_when(st == "VT"~7,
                      st == "NC"~6.125,
                      st == "CO"~7.45,
                      st == "OR"~8,
                      T~8) / 1.25
  
  p <- tolower(paste0(gsub(" ", "-", spec), "_", st, "_obs"))
  
  if (!file.exists(file.path("plots", paste0(p, ".svg"))) | refresh.obs.plts) {
    p <- get.pres.abs.plt(data.sf, states.sf, st, spec, pres="Presence", 
                           title=paste0(spec, " Observations in ", st)) %>%
      plot.to.svg(p=., st, 
                  plt.name=p, 
                  wh=list(w=width, h=height)) 
  } else {
    p <- file.path("plots", paste0(p, ".svg"))
  }
  plt.container <- div(
    id=tolower(paste0(st, "_", gsub(" ", "-", spec), "_obs_plots")),
    style=ifelse(i==1, "", "display:none;"),
    div(
      style="margin:5px;",
      tags$img(src=p)
    )
  )
  
})


obs.plt.selections <- htmltools::div(
  htmltools::tags$script(
    '$(document).ready(function(){
        $("#obs_state_selector").change(function(){
          var selectedState = $(this).val();
          var selectedSpecies = $("#obs_species_selector").val();
          // Hide all raster plots
          $("[id$=_obs_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_obs_plots").show();
        });
        $("#obs_species_selector").change(function(){
          var selectedState = $("#obs_state_selector").val();
          var selectedSpecies = $(this).val();
          // Hide all raster plots
          $("[id$=_obs_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_obs_plots").show();
        });
      });'
  ),
  htmltools::tags$div(
    style="display:flex; flex-direction:row;",
    htmltools::tags$div(
      style="margin-right:5px;",
      htmltools::tags$select(id='obs_state_selector',
                             style="font-size:17px;",
                             lapply(states, function(s) {
                               htmltools::tags$option(value=tolower(s), s)
                             }))
    ),
    htmltools::tags$div(
      htmltools::tags$select(id='obs_species_selector',
                             style="font-size:17px;",
                             lapply(species, function(s) {
                               htmltools::tags$option(value=tolower(gsub(" ", "-", s)), s)
                             }))
    )
  ),
  obs.list
)

sample.obs.img <- (get.pres.abs.plt(data.sf, states.sf, "OR", "Wild Turkey", 
                                  pres="Presence", 
                                  title=NULL) +
  theme(plot.margin = margin(r=80, l=15, t=0, b=0, unit="pt"))) %>%
  plot.to.svg(p=., "OR", 
              plt.name="obs_sample",
              wh=list(h=4, w=5.5))



# MODEL COVARIATES --------------------------------------------------------

# Load prepared explanatory rasters by state

covariate.list <- states %>%
  set_names %>%
  purrr::map(~{
    f <- paste0("../data/final_rasters/", .x, ".tif")
    rast(f) %>%
      project(., crs(pred.rasters[[1]][[1]][[1]]))
  })

covariate.list$OR$Summer_NDVI_OR

sample.cov.img <- plot.raster.gg(covariate.list$OR$Summer_NDVI_OR, 
                                "Summer_NDVI_OR", scale.c=T, 
                                title="",
                                legend.title="Summer NDVI") %>%
  plot.to.svg(p=., "OR", 
              plt.name="covariate_sample",
              wh=list(h=4, w=5.5))

# PSEUDO-ABSENCE PLOTS ----------------------------------------------------
refresh.pa.plts <- F
pa.list <- purrr::map(1:nrow(spec.state), function(i) {
  st <- spec.state[i,]$state
  spec <- spec.state[i,]$species
  r <- bioclim.list[[i]]$bioclim %>%
    project(., crs(pred.rasters[[1]][[1]][[1]]))
  qnt.5 <- bioclim.list[[i]]$qnt.5
  fpar <- bioclim.list[[i]]$final.pseudo.absence.region %>%
    project(., crs(pred.rasters[[1]][[1]][[1]]))
  width <- case_when(st == "VT"~7.5,
                     st == "NC"~14,
                     st == "CO"~10.5,
                     st == "OR"~10.6,
                     T~8) / 1.65
  height <- case_when(st == "VT"~7,
                      st == "NC"~6.125,
                      st == "CO"~7.45,
                      st == "OR"~8,
                      T~8) / 1.65
  p1 <- tolower(paste0(gsub(" ", "-", spec), "_", st, "_bioclim"))
  p2 <- tolower(paste0(gsub(" ", "-", spec), "_", st, "_pa1"))
  p3 <- tolower(paste0(gsub(" ", "-", spec), "_", st, "_pa2"))
  p4 <- tolower(paste0(gsub(" ", "-", spec), "_", st, "_pa_points"))
  
  if (!file.exists(file.path("plots", paste0(p1, ".svg"))) | refresh.pa.plts) {
    p1 <- plot.raster.gg(r, "sum", legend.title="BIOCLIM Sum", scale.c=T,
                         title=paste0("BIOCLIM Suitability for\n", 
                                      spec, " in ", st)) %>%
      plot.to.svg(p=., st, 
                  plt.name=p1, 
                  wh=list(w=width, h=height)) 
  } else {
    p1 <- file.path("plots", paste0(p1, ".svg"))
  }
  if (!file.exists(file.path("plots", paste0(p2, ".svg"))) | refresh.pa.plts) {
    p2 <- plot.raster.gg(r < qnt.5, "sum", scale.c=F,
                         title = paste0("BIOCLIM Suitable sample regions\nfor", 
                                        spec, " in ", st),
                         legend.title="BIOCLIM Suitable\nPseudo-absence\nRegions") %>%
      plot.to.svg(p=., st, 
                  plt.name=p2, 
                  wh=list(w=width, h=height)) 
  } else {
    p2 <- file.path("plots", paste0(p2, ".svg"))
  }
  if (!file.exists(file.path("plots", paste0(p3, ".svg"))) | refresh.pa.plts) {
    p3 <- plot.raster.gg(fpar == 1, "sum", scale.c=F,
                         title =paste0("Final Suitable sample regions\nfor", 
                                       spec, " in ", st),
                         legend.title="Suitable Sampling\nRegions") %>%
      plot.to.svg(p=., st, 
                  plt.name=p3, 
                  wh=list(w=width, h=height)) 
  } else {
    p3 <- file.path("plots", paste0(p3, ".svg"))
  }
  
  if (!file.exists(file.path("plots", paste0(p4, ".svg"))) | refresh.pa.plts) {
    p4 <- get.pres.abs.plt(data.sf, states.sf, st, spec, pres="both", 
                     title=paste0("Observations and Pseudo-Absence Points\n",
                                  "for the ", spec, " in ", st)) %>%
      plot.to.svg(p=., st, 
                  plt.name=p4, 
                  wh=list(w=width, h=height)) 
  } else {
    p4 <- file.path("plots", paste0(p4, ".svg"))
  }
  
  plots <- div(
    id=tolower(paste0(st, "_", gsub(" ", "-", spec), "_bioclim_plots")),
    style=paste0(ifelse(i==1, "display:flex", "display:none;"), "flex-direction:column;"),
    div(
      style="display:flex; flex-direction:row;",
      tags$div(
        style="margin:5px;",
        tags$img(src=p1)
      ),
      tags$div(
        style="margin:5px;",
        tags$img(src=p2)
      )
    ),
    div(
      style="display:flex; flex-direction:row;",
      tags$div(
        style="margin:5px;",
        tags$img(src=p3)
      ),
      tags$div(
        style="margin:5px;",
        tags$img(src=p4)
      )
    )
  )
})

pa.plt.selections <- htmltools::div(
  htmltools::tags$script(
    '$(document).ready(function(){
        $("#bioclim_state_selector").change(function(){
          var selectedState = $(this).val();
          var selectedSpecies = $("#bioclim_species_selector").val();
          // Hide all raster plots
          $("[id$=_bioclim_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_bioclim_plots").show();
        });
        $("#bioclim_species_selector").change(function(){
          var selectedState = $("#bioclim_state_selector").val();
          var selectedSpecies = $(this).val();
          // Hide all raster plots
          $("[id$=_bioclim_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_bioclim_plots").show();
        });
      });'
  ),
  htmltools::tags$div(
    style="display:flex; flex-direction:row;",
    htmltools::tags$div(
      style="margin-right:5px;",
      htmltools::tags$select(id='bioclim_state_selector',
                             style="font-size:17px;",
                             lapply(states, function(s) {
                               htmltools::tags$option(value=tolower(s), s)
                             }))
    ),
    htmltools::tags$div(
      htmltools::tags$select(id='bioclim_species_selector',
                             style="font-size:17px;",
                             lapply(species, function(s) {
                               htmltools::tags$option(value=tolower(gsub(" ", "-", s)), s)
                             }))
    )
  ),
  pa.list
)

sample.pa.img <- bioclim.list[[24]]$final.pseudo.absence.region %>%
  project(., crs(pred.rasters[[1]][[1]][[1]]))
sample.pa.img <- plot.raster.gg(sample.pa.img == 1, "sum", scale.c=F, title="",
                 legend.title="Suitable\nSampling\nRegions") %>%
  plot.to.svg(p=., "OR", 
              plt.name="pa_sample",
              wh=list(h=4, w=5.5))

# PREDICTION PLOTS --------------------------------------------------------

# Create plots for each combination
pred.plts.svg <- get.object(
  purrr::map_df(1:nrow(spec.state), function(i) {
    st <- spec.state[i,]$state
    spec <- spec.state[i,]$species
    model.types <- spec.state.mods[species == spec & state == st]$model.type
    cat(st, spec, "\n")
    r <- model.types %>% 
      set_names() %>% 
      purrr::map(~pred.rasters[[.x]][[spec]][[st]])
    plots <- model.types %>%
      set_names %>%
      purrr::map(~{
        mt.case <- get.model.type.formatted(.x)
        plot.raster.gg(
          r=r[[.x]],
          fill=names(r[[.x]]),
          scale.c=T,
          legend.title=case_when(names(r[[.x]]) == "trend" ~ "trend",
                                 names(r[[.x]]) == "prob" ~ "probability",
                                 T ~ "prediction"),
          title=glue::glue("{mt.case} prediction for \n{spec}s in {st}")
        )
      })
    p <- ggpubr::ggarrange(plotlist=plots, ncol=4, nrow=2)
    p.svg <- plot.to.svg(p, st, plt.name=tolower(paste0(gsub(" ", "-", spec), "_", st)))
    data.table(species=spec, state=st, plot.svg=p.svg)
  }),
  file.name="pred_plots_svg.rds",
  obj.path=getwd()
)


pred.plt.selections <- htmltools::div(
  htmltools::tags$script(
    '$(document).ready(function(){
        $("#pred_state_selector").change(function(){
          var selectedState = $(this).val();
          var selectedSpecies = $("#pred_species_selector").val();
          // console.log(selectedState + " & " + selectedSpecies);
          // Hide all raster plots
          $("[id$=_pred_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_pred_plots").show();
        });
        $("#pred_species_selector").change(function(){
          var selectedState = $("#pred_state_selector").val();
          var selectedSpecies = $(this).val();
          // console.log(selectedState + " & " + selectedSpecies);
          // Hide all raster plots
          $("[id$=_pred_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_pred_plots").show();
        });
      });'
  ),
  htmltools::tags$div(
    style="display:flex; flex-direction:row;",
    htmltools::tags$div(
      style="margin-right:5px;",
      htmltools::tags$select(id='pred_state_selector',
                             style="font-size:17px;",
                             lapply(states, function(s) {
                               htmltools::tags$option(value=s, s)
                             }))
    ),
    htmltools::tags$div(
      style="margin-right:5px;",
      htmltools::tags$select(id='pred_species_selector',
                             style="font-size:17px;",
                             lapply(species, function(s) {
                               htmltools::tags$option(value=gsub(" ", "-", s), s)
                             }))
    )
  ),
  purrr::map(1:nrow(spec.state), function(i) {
    st <- spec.state.mods[i,]$state
    spec <- spec.state.mods[i,]$species
    d <- pred.plts.svg[state == st & species == spec]
    p.svg <- d$plot.svg %>% basename()
    # img.src <- paste0("https://raw.githubusercontent.com/benton-tripp/",
    #                   "presence-only-sdm/main/app/plots/", p.svg)
    img.src <- paste0("plots/", p.svg)
    htmltools::tags$div(
      id=paste0(st, "_", gsub(" ", "-", spec), "_pred_plots"),
      style=paste0("padding:5px; overflow:auto; width:100%;", 
                   "min-height:1400px; ", 
                   ifelse(i==1, "", "display:none;")),
      htmltools::tags$div(
        style="margin:5px;",
        htmltools::tags$img(
          src=img.src, 
          style="max-height:700px; width=auto; height=auto;"
        )
      )
    )
  })
)

sample.pred.img <- plot.raster.gg(
  r=pred.rasters[["maxent"]][["Wild Turkey"]][["OR"]],
  fill="prob",
  scale.c=T,
  legend.title="Probability",
  title="") %>%
  plot.to.svg(., "OR", plt.name="prediction_sample",
              wh=list(h=4, w=5.5))


# CREATE UI ---------------------------------------------------------------


ui.head <- htmltools::tagList(
  tags$meta(
    `http-equiv`="Content-Type",
    content="text/html; charset=utf-8"
  ),
  tags$script(type="application/shiny-singletons"),
  tags$script(
    type="application/html-dependencies",
    paste0(
      "jquery[3.6.0];",
      "shiny-css[1.7.4];",
      "font-awesome[6.4.0];",
      "bootstrap[3.4.1]"
    )
  ),
  tags$script(src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.0/jquery.min.js"),
  tags$script(
    src="https://cdn.jsdelivr.net/npm/bootstrap@3.4.1/dist/js/bootstrap.min.js",
    rel="stylesheet",
    integrity="sha384-aJ21OjlMXNL5UyIl/XNwTMqvzeRMZH2w8c5cRVpzpU8Y5bApTppSuUkhZXN0VxHd",
    crossorigin="anonymous"
  ),
  tags$link(
    href="shiny.min.css", 
    rel="stylesheet"
  ),
  tags$link(
    href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0/css/all.min.css",
            rel="stylesheet"
  ),
  tags$link(
    href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.4.0//css/v4-shims.min.css", 
    rel="stylesheet"
  ),
  tags$meta(name="viewport", content="width=device-width, initial-scale=1"),
  tags$link(
    rel="stylesheet", 
    href="https://cdn.jsdelivr.net/npm/bootstrap@3.4.1/dist/css/bootstrap.min.css",
    integrity="sha384-HSMxcRTRxnN+Bdg0JdbxYKrThecOKuH5zCYotlSAcp1+c8xmyTe9GYg1l9a69psu",
    crossorigin="anonymous"
  ),
  tags$script(
    rel="stylesheet",
    href="https://cdn.jsdelivr.net/npm/bootstrap@3.4.1/dist/css/bootstrap-theme.min.css",
    integrity="sha384-6pzBo3FDv/PJ8r2KRkGHifhEocL+1X2rVCTTkUfGk7/0pbek5mMa1upzvWbrUbOZ",
    crossorigin="anonymous"
  ),
  tags$title("SDM Benchmark Study Results"),
  # tags$link(rel="stylesheet", type="text/css", href="bootstrap.css"),
  tags$link(rel="stylesheet", type="text/css", href="styles.css"),
  tags$script(src="scripts.js")
  # tags$script(src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.4/jquery.min.js")
)

ui.nav <- htmltools::tags$nav(
  class="navbar navbar-default navbar-fixed-top",
  role="navigation",
  div(
    class="container-fluid",
    div(
      class="navbar-header active shared-menu-item",
      tags$a(
        href="#tab-5466-1",
        `data-toggle`="tab",
        `data-bs-toggle`="tab",
        `data-value`="Home",
        tags$span(
          class="navbar-brand",
          "SDM Benchmark Study Results"
        )
      )
    ),
    tags$ul(
      class="nav navbar-nav navbar-right",
      id="mainPage",
      `data-tabsetid`="1639",
      tags$li(
        class="dropdown",
        tags$a(
          id="dropDownMenu",
          href="#",
          type="button",
          class="btn btn-default dropdown-toggle",
          `data-toggle`="dropdown",
          `data-value`="mainMenu",
          `aria-haspopup`="true",
          `aria-expanded`="false",
          tags$i(
            `aria-label`="bars icon",
            class="fas fa-bars fa-fw",
            role="presentation"
          ),
          tags$span("Menu"),
          tags$b(
            class="caret"
          )
        ),
        tags$ul(
          class="dropdown-menu",
          `aria-labelledby`="dropDownMenu",
          `data-tabsetid`="5466",
          tags$li(
            class="active shared-menu-item",
            tags$a(
              class="dropdown-item",
              href="#tab-5466-1",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Home",
              tags$b("Home")
            )
          ),
          tags$li(
            class="shared-menu-item",
            tags$a(
              class="dropdown-item",
              href="#tab-5466-2",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Species Observations",
              "Species Observations"
            )
          ),
          tags$li(
            class="shared-menu-item",
            tags$a(
              class="dropdown-item",
              href="#tab-5466-3",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Model Covariates",
              "Model Covariates"
            )
          ),
          tags$li(
            class="shared-menu-item",
            tags$a(
              class="dropdown-item",
              href="#tab-5466-4",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Pseudo-Absence Selection",
              "Pseudo-Absence Selection"
            )
          ),
          tags$li(
            class="shared-menu-item",
            tags$a(
              class="dropdown-item",
              href="#tab-5466-5",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Species Predictions",
              "Species Predictions"
            )
          ),
          hr(),
          tags$li(
            tags$a(
              class="dropdown-item",
              href="https://benton-tripp.github.io/",
              tags$span("Benton's Portfolio/Blog")
            )
          )
        )
      )
    )
  )
)

ui.main <- div(
  class = "container-fluid",
  div(class = "row"),
  div(
    class = "tab-content", 
    `data-tabsetid` = "1639",
    div(
      class = "tab-pane active", 
      `data-value` = "Home", 
      id = "tab-5466-1",
      div(
        class = "main-area-container",
        h1("Species Distribution Benchmark Study Results"),
        h4("Benton Tripp", style="color:#555555;"),
        div(
          id="homePageMain",
          div(
            id="homePageContainer",
            tags$ul(
              class = "nav nav-tabs",
              role="tablist",
              style="list-style-type:none; display:flex; flex-direction:column;",
              tags$li(
                class="shared-menu-item",
                style="margin:5px;",
                tags$a(
                  class="link-header-parent",
                  href = "#tab-5466-2",
                  `data-toggle`="tab",
                  tags$span(
                    class="link-header",
                    "Species Observations"
                  ),
                  tags$img(
                    style="width:250px; height:auto;",
                    src=sample.obs.img
                  )
                )
              ),
              tags$li(
                class="shared-menu-item",
                style="margin:5px;",
                tags$a(
                  class="link-header-parent",
                  href = "#tab-5466-3",
                  `data-toggle` = "tab",
                  tags$span(
                    class="link-header",
                    "Model Covariates"
                  ),
                  tags$img(
                    style="width:250px; height:auto;",
                    src=sample.cov.img 
                  )
                )
              ),
              tags$li(
                class="shared-menu-item",
                style="margin:5px;",
                tags$a(
                  class="link-header-parent",
                  href = "#tab-5466-4",
                  `data-toggle` = "tab",
                  tags$span(
                    class="link-header",
                    "Pseudo-Absence Selection"
                  ),
                  tags$img(
                    style="width:250px; height:auto;",
                    src=sample.pa.img 
                  )
                )
              ),
              tags$li(
                class="shared-menu-item",
                style="margin:5px;",
                tags$a(
                  class="link-header-parent",
                  href="#tab-5466-5",
                  `data-toggle`="tab",
                  `data-value`="Species Predictions",
                  tags$span(
                    class="link-header",
                    "Species Predictions"
                  ),
                  tags$img(
                    src=sample.pred.img, 
                    style="width:250px; height:auto;",
                  )
                )
              )
            )
          ),
          div(
            id="mainReportContainer",
            h3("Project Overview and Results"),
            tags$ul(
              style="list-style-type:none; display:flex; flex-direction:column;",
              tags$li(
                tags$a(
                  href="https://raw.githubusercontent.com/benton-tripp/benton-tripp.github.io/main/_docs/lit_review.pdf",
                  target="_blank",
                  HTML('<i class="fa-solid fa-book" style="margin-right:3px;"></i>'),
                  tags$span("Literature Review"),
                )
              ),
              tags$li(
                tags$a(
                  href="https://raw.githubusercontent.com/benton-tripp/presence-only-sdm/main/docs/FinalSummary.pdf",
                  target="_blank",
                  HTML('<i class="fa-solid fa-file" style="margin-right:3px;"></i>'),
                  tags$span("Final Results and Summary"),
                )
              ),
              tags$li(
                tags$a(
                  href="https://github.com/benton-tripp/presence-only-sdm",
                  target="_blank",
                  HTML('<i class="fa-brands fa-github" style="margin-right:3px;"></i>'),
                  tags$span("Github Repository")
                )
              ),
              tags$li(
                tags$p(
                  "<Overview here>"
                )
              )
            ),
            h3("Project Walk-Through"),
            tags$ol(
              style="display:flex; flex-direction:column;",
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/sdm-benchmark-study-part-1-data-preparation.html",
                  target="_blank",
                  tags$span("Data Preparation")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-09-17-sdm-benchmark-study-part-2-exploratory-analysis.html",
                  target="_blank",
                  tags$span("Exploratory Analysis")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-09-27-sdm-benchmark-study-part-3-more-preprocessing-and-eda.html",
                  target="_blank",
                  tags$span("More Pre-Processing and EDA")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-10-05-sdm-benchmark-study-part-4-baseline-species-distribution-models.html",
                  target="_blank",
                  tags$span("Baseline Species Distribution Models")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-10-13-sdm-benchmark-study-part-5-fitting-and-testing-ipp-models.html",
                  target="_blank",
                  tags$span("Fitting and Testing Inhomogeneous Poisson Process Models ")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-10-18-sdm-benchmark-study-part-6-resampling-pseudoabsence-points.html",
                  target="_blank",
                  tags$span("Re-Sampling Pseudo-Absence Points Using Iterative Sampling and BIOCLIM ")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-10-26-sdm-benchmark-study-part-7-ipp-models-with-updated-data.html",
                  target="_blank",
                  tags$span("Fitting and Testing Inhomogeneous Poisson Process Models with Updated Data ")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-11-03-sdm-benchmark-study-part-8-fitting-and-testing-maxent-models.html",
                  target="_blank",
                  tags$span("Fitting and Testing MaxEnt Models ")
                )
              ),
              tags$li(
                tags$a(
                  href="https://benton-tripp.github.io/posts/2023-11-23-sdm-benchmark-study-part-9-fitting-and-testing-ml-models.html",
                  target="_blank",
                  tags$span("Fitting and Testing Common ML Models ")
                )
              )
            )
          )
        )
      )
    ),
    div(
      class = "tab-pane", 
      `data-value` = "Species Observations", 
      id = "tab-5466-2",
      div(
        class = "main-area-container",
        h1("Species Observations"),
        obs.plt.selections
      )
    ),
    div(
      class = "tab-pane", 
      `data-value` = "Model Covariates", 
      id = "tab-5466-3",
      div(
        class = "main-area-container",
        h1("Model Covariates"),
        tags$i("*Prior to pre-processing and feature engineering")
      )
    ),
    div(
      class = "tab-pane", 
      `data-value` = "Pseudo-Absence Selection", 
      id = "tab-5466-4",
      div(
        class = "main-area-container",
        HTML(paste0("<h1 style='margin-bottom:10px;'>BIOCLIM Suitability, ",
                    "Pseudo-Absence Sample Regions,",
                    "<br>and Pseudo-Absence Selections")),
        pa.plt.selections
      )
    ),
    div(
      class = "tab-pane", 
      `data-value` = "Species Predictions", 
      id = "tab-5466-5",
      div(
        class = "main-area-container",
        h1("Species Predictions"),
        pred.plt.selections
      )
    )
  )
)


html.output <- tags$html(
  htmltools::HTML(glue::glue("<head>{ui.head}</head>")),
  tags$body(
    ui.nav,
    ui.main
  )
)
readr::write_file(glue::glue("{html.output}"), "SDM-Benchmark-Study-Results.html")
