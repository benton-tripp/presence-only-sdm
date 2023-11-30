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
    tabPanel(
      "Model Covariates",
      div(
        class="main-area-container",
        h1("Model Covariates"),
        tags$i("*Prior to pre-processing and feature engineering")
      )
    ),
    tabPanel(
      "Pseudo-Absence Selection",
      div(
        class="main-area-container",
        h1("Pseudo-Absence Selection"),
        div(
          class="main-area-container",
          tabsetPanel(
            id="paTabsetPanel",
            tabPanel(
              "BIOCLIM",
              div(
                h2("BIOCLIM Suitability")
              )
            ),
            tabPanel(
              "Sample Regions",
              h2("Final Sample Regions")
            ),
            tabPanel(
              "Pseudo-Absence Points",
              h2("Pseudo-Absence Points")
            )
          )
        )
      )
    ),
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
      class="navbar-header",
      tags$span(
        class="navbar-brand",
        "SDM Benchmark Study Results"
      )
    ),
    tags$ul(
      class="nav navbar-nav shiny-tab-input",
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
          "Menu",
          tags$b(
            class="caret"
          )
        ),
        tags$ul(
          class="dropdown-menu",
          `aria-labelledby`="dropDownMenu",
          `data-tabsetid`="5466",
          tags$li(
            class="active",
            tags$a(
              class="dropdown-item",
              href="#tab-5466-1",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Species Observations",
              "Species Observations"
            )
          ),
          tags$li(
            tags$a(
              class="dropdown-item",
              href="#tab-5466-2",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Model Covariates",
              "Model Covariates"
            )
          ),
          tags$li(
            tags$a(
              class="dropdown-item",
              href="#tab-5466-3",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Pseudo-Absence Selection",
              "Pseudo-Absence Selection"
            )
          ),
          tags$li(
            tags$a(
              href="#tab-5466-4",
              `data-toggle`="tab",
              `data-bs-toggle`="tab",
              `data-value`="Species Predictions",
              "Species Predictions"
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
      `data-value` = "Species Observations", 
      id = "tab-5466-1",
      div(
        class = "main-area-container",
        h1("Species Observations")
      )
    ),
    div(
      class = "tab-pane", 
      `data-value` = "Model Covariates", 
      id = "tab-5466-2",
      div(
        class = "main-area-container",
        h1("Model Covariates"),
        tags$i("*Prior to pre-processing and feature engineering")
      )
    ),
    div(
      class = "tab-pane", 
      `data-value` = "Pseudo-Absence Selection", 
      id = "tab-5466-3",
      div(class = "main-area-container",
          h1("Pseudo-Absence Selection"),
          div(class = "tabbable",
              tags$ul(
                class = "nav nav-tabs shiny-tab-input", 
                id = "paTabsetPanel", 
                `data-tabsetid` = "5409",
                tags$li(
                  class = "active",
                  tags$a(
                    href = "#tab-5409-1", 
                    `data-toggle` = "tab", 
                    `data-bs-toggle` = "tab", 
                    `data-value` = "BIOCLIM", "BIOCLIM")
                ),
                tags$li(
                  tags$a(
                    href = "#tab-5409-2", 
                    `data-toggle` = "tab", 
                    `data-bs-toggle` = "tab", 
                    `data-value` = "Sample Regions", 
                    "Sample Regions")
                ),
                tags$li(
                  tags$a(
                    href = "#tab-5409-3", 
                    `data-toggle` = "tab", 
                    `data-bs-toggle` = "tab",
                    `data-value` = "Pseudo-Absence Points", 
                    "Pseudo-Absence Points")
                )
              ),
              div(
                class = "tab-content", 
                `data-tabsetid` = "5409",
                div(
                  class = "tab-pane active", 
                  `data-value` = "BIOCLIM", 
                  id = "tab-5409-1",
                  div(
                    h2("BIOCLIM Suitability")
                  )
                ),
                div(
                  class = "tab-pane", 
                  `data-value` = "Sample Regions", 
                  id = "tab-5409-2",
                  h2("Final Sample Regions")
                ),
                div(
                  class = "tab-pane", 
                  `data-value` = "Pseudo-Absence Points", 
                  id = "tab-5409-3",
                  h2("Pseudo-Absence Points")
                )
              )
          )
      )
    ),
    div(
      class = "tab-pane", 
      `data-value` = "Species Predictions", 
      id = "tab-5466-4",
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
