library(shiny)
library(shinyIncubator)
library(RMySQL)
options(stringsAsFactors=F)

shinyUI(
  fluidPage(theme="css/united.bootstrap.min.css",
            title="CaVaDa",
            progressInit(),
            tagList(
              singleton(
                tags$head(
                  tags$script(src = "js/jsCodeMessage.js"), 
                  #tags$script(src = "js/disable_back.js"), 
                  tags$style(type="text/css", '.help-block {color: #A0A0A0; margin-bottom:5px;}'),
                  tags$style(type="text/css", 'label[for=offsetLimit] {margin-left:20px; display:inline-block;}'),
                  tags$style(type="text/css", '#offsetLimit {width:50px; margin-bottom:0px; margin-left:5px; margin-right:5px;}'),
                  tags$style(type="text/css", '.thin {min-height:20px !important;}'),
                  tags$style(type="text/css",".highlight {background-color: rgba(232, 101, 55, 0.6) !important}"),
                  tags$style(type="text/css",'.row-fluid [class*="span"] {min-height: 20px !important}'),
                  tags$style(type="text/css",'.details {background-color: #D1CFD0 !important;border: 1px solid #A19B9E;}'),
                  HTML('<style type="text/css"> td > strong {margin-left:5px;}</style>')
                )
              )
            ),
            uiOutput('uiControl'),
            uiOutput('uiContents')
  )
)
