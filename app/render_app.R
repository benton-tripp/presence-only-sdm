setwd("C:/Users/bento/gis630/app")
load(file="../artifacts/summary_img.RData")
library(data.table)
library(dplyr)
library(purrr)
library(terra)
library(caret)
library(sf)
library(spatstat)
library(shiny)

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
  width <- case_when(st == "VT"~5.25*4/2,
                     st == "NC"~16*4/2,
                     st == "CO"~12*4/2,
                     st == "OR"~8*4/2,
                     T~8*4/2)
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

plot.to.svg <- function(p, st, encode=F, plt.dir="plots", plt.name="plt") {
  # Create a temporary file
  tmpfile <- tempfile(fileext = ".svg")
  # Save the ggplot to this file
  wh <- plt.width.height(st)
  ggsave(filename = tmpfile, p, width=wh$w, height=wh$h)
  if (encode) {
    # Use img() function from htmltools to display it within a div
    img.src <- paste0("data:image/svg+xml;base64,", base64enc::base64encode(tmpfile))
  } else {
    img.src <- file.path(plt.dir, paste0(plt.name, ".svg"))
    file.copy(tmpfile, img.src)
  }
  file.remove(tmpfile)
  img.src
}

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
          console.log(selectedState + " & " + selectedSpecies);
          // Hide all raster plots
          $("[id$=_pred_plots]").hide();
          // Show the selected raster plot
          $("#" + selectedState + "_" + selectedSpecies + "_pred_plots").show();
        });
        $("#pred_species_selector").change(function(){
          var selectedState = $("#pred_state_selector").val();
          var selectedSpecies = $(this).val();
          console.log(selectedState + " & " + selectedSpecies);
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
    img.src <- paste0("https://raw.githubusercontent.com/benton-tripp/",
                      "presence-only-sdm/main/app/plots/", p.svg)
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

ui <- navbarPage(
  title="SDM Benchmark Study Results",
  id="mainPage",
  position="fixed-top",
  theme="boostrap.css",
  header=htmltools::tags$head(
    includeScript("scripts.js"),
    includeCSS("styles.css"),
    htmltools::tags$script(
      src="https://ajax.googleapis.com/ajax/libs/jquery/3.6.4/jquery.min.js"
    )
  ),
  selected="Species Observations",
  navbarMenu(
    title="Menu",
    menuName="mainMenu",
    icon=icon("bars", verify_fa=F),
    tabPanel(
      "Species Observations",
      div(
        class="main-area-container",
        h1("Species Observations")
      )
    ),
    tabPanel("Covariates"),
    tabPanel("Pseudo-Absence Selection"),
    tabPanel(
      "Species Predictions",
      div(
        class="main-area-container",
        h1("Species Predictions"),
        pred.plt.selections
      )
    )
  )
)

server <- function(input, output) {
  
}

shiny::shinyApp(ui=ui, server=server)
